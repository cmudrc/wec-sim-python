"""Gate the corrected OSWEC scene from independently advanced Python motion."""

import os
from pathlib import Path
from xml.etree import ElementTree
from zipfile import ZipFile

import numpy as np
import pytest
from scipy.optimize import brentq
from scipy.spatial import cKDTree
import trimesh

from wecsim import RegularWave, WEC, WorldPoint
from wecsim.environment.waveClass import WaveClass
from wecsim.hydrodynamics.nonlinearHydro import regular_wave_mesh_pressures
from wecsim.visualization.paraviewClass import ParaviewClass


ROOT = Path(__file__).parent
FIXTURE = ROOT / "fixtures/matlab/oswec_nonlinear_viz_corrected"
HYDRO = ROOT / "fixtures/matlab/oswec_multi_wave/oswec.h5"
FLAP = ROOT / "test_objects/test_bodyclass/testData/body_3_test/flap.stl"
BASE = ROOT / "test_objects/test_bodyclass/testData/body_3_test/base.stl"


def _maximum(actual, reference, limit, label):
    gap = float(np.max(np.abs(np.asarray(actual) - np.asarray(reference))))
    assert gap < limit, f"{label}: {gap:.9g} >= {limit:.9g}"


def _vtk(path_or_xml):
    root = (ElementTree.fromstring(path_or_xml)
            if isinstance(path_or_xml, bytes)
            else ElementTree.parse(path_or_xml).getroot())
    piece, = root.findall("./PolyData/Piece")
    points = np.fromstring(piece.findtext("./Points/DataArray"), sep=" ").reshape(-1, 3)
    polys = piece.find("Polys")
    cells = np.fromstring(polys.find("./DataArray[@Name='connectivity']").text,
                          sep=" ", dtype=int)
    offsets = np.fromstring(polys.find("./DataArray[@Name='offsets']").text,
                            sep=" ", dtype=int)
    fields = {item.get("Name"): np.fromstring(item.text, sep=" ")
              for item in piece.findall("./CellData/DataArray")}
    return points, cells, offsets, fields


def _run_independent_case():
    wec = WEC("Corrected nonlinear OSWEC validation")
    flap = wec.body("flap", HYDRO, mass=127000, inertia=(1.85e6,) * 3,
                    geometry_file=FLAP, nonlinear_hydro="instantaneous")
    base = wec.fixed_body("base", center_gravity=(0, 0, -10.9),
                          mass=999, inertia=(1, 1, 1))
    wec.fixed_hinge(flap, base, location=WorldPoint(0, 0, -10),
                    pto_location=WorldPoint(0, 0, -8.9))
    return wec.run(RegularWave(2.5, 8), dt=0.1, end_time=120, ramp_time=40)


def _compare_scene_frame(frame, wave_path, flap_path, base_path, source_file):
    actual = _vtk(wave_path)
    expected = _vtk(source_file(f"waves/waves_{frame}.vtp"))
    _maximum(actual[0], expected[0], 1.01e-5, "wave vertices")
    np.testing.assert_array_equal(actual[1], expected[1])
    np.testing.assert_array_equal(actual[2], expected[2])
    for folder, stem, path, vertex_limit in (
        ("body1_flap", "flap", flap_path, .0025),
        ("body2_base", "base", base_path, 1e-5),
    ):
        actual = _vtk(path)
        expected = _vtk(source_file(f"{folder}/{stem}_{frame}.vtp"))
        actual_points, actual_cells, actual_offsets, actual_fields = actual
        source_points, source_cells, source_offsets, source_fields = expected
        np.testing.assert_array_equal(actual_offsets, source_offsets)
        faces = actual_points[actual_cells.reshape(-1, 3)]
        source_faces = source_points[source_cells.reshape(-1, 3)]
        distance, match = cKDTree(faces.mean(axis=1)).query(
            source_faces.mean(axis=1))
        assert len(np.unique(match)) == len(match)
        assert np.max(distance) < vertex_limit
        vertex_error = np.linalg.norm(
            source_faces[:, :, None, :] - faces[match][:, None, :, :],
            axis=3,
        )
        assert np.max(np.min(vertex_error, axis=2)) < vertex_limit
        assert actual_fields.keys() == source_fields.keys()
        _maximum(actual_fields["CellArea"][match],
                 source_fields["CellArea"], 1e-3, f"{stem} cell areas")
        if stem == "flap":
            combined = (actual_fields["HydrostaticPressure"][match]
                        + actual_fields["WavePressureNonLinear"][match])
            source_combined = (source_fields["HydrostaticPressure"]
                               + source_fields["WavePressureNonLinear"])
            _maximum(combined, source_combined, 13,
                     "VTP combined pressure")
            _maximum(actual_fields["WavePressureLinear"][match],
                     source_fields["WavePressureLinear"], 3,
                     "VTP linear pressure")


def test_corrected_oswec_full_motion_pressure_and_scene(tmp_path):
    with np.load(FIXTURE / "history.npz", allow_pickle=False) as archive:
        source = {key: archive[key] for key in archive.files}
    time = source["time"]
    np.testing.assert_allclose(time, np.arange(1201) * 0.1, rtol=0, atol=1e-11)
    np.testing.assert_array_equal(source["source_dt"], [0.0125])
    np.testing.assert_array_equal(source["stability_dt"], [0.00625])

    # A second independently run MATLAB arm changes only the numerical step.
    for field, limits in (("pose", {0: .0007, 2: .0003, 4: .00015}),
                          ("velocity", {0: .0006, 2: .0003, 4: .00012})):
        for axis, limit in limits.items():
            _maximum(source[field][:, axis], source[f"stable_{field}"][:, axis],
                     limit, f"source stability {field} axis {axis}")

    result = _run_independent_case()
    np.testing.assert_allclose(result.time, time, rtol=0, atol=1e-11)
    pose = result.bodies["flap"].position
    speed = result.bodies["flap"].velocity
    for field, actual, limits in (
        ("pose", pose, {0: .001, 2: .0005, 4: .00025}),
        ("velocity", speed, {0: .001, 2: .0005, 4: .0002}),
    ):
        for axis in range(6):
            _maximum(actual[:, axis], source[field][:, axis],
                     limits.get(axis, 1e-10), f"flap {field} axis {axis}")
    for name, field in (("position", "basePose"),
                        ("velocity", "baseVelocity")):
        _maximum(getattr(result.bodies["base"], name), source[field],
                 1e-10, f"base {name}")
    pto = result.ptos["hinge"]
    for label, actual, field, limit in (
        ("stroke", pto.stroke, "ptoPosition", .00025),
        ("speed", pto.velocity, "ptoVelocity", .0002),
        ("force", pto.force, "ptoForce", 1e-5),
        ("power", pto.absorbed_power, "ptoPower", 1e-6),
    ):
        _maximum(actual, source[field][:, 4], limit, f"PTO {label}")

    depth = 10.9
    gravity = 9.81
    omega = 2 * np.pi / 8
    wave_number = brentq(
        lambda k: gravity * k * np.tanh(k * depth) - omega**2,
        1e-12, 1,
    )
    _maximum(wave_number, source["wave_number"], 1e-12, "wave number")
    np.testing.assert_array_equal(source["water_depth"], [depth])
    mesh = trimesh.load_mesh(FLAP, process=True)
    assert len(mesh.faces) == len(source["areas"]) == 1042
    _maximum(mesh.area_faces, source["areas"], 1e-9, "flap facet areas")
    pressure = regular_wave_mesh_pressures(
        mesh.vertices, mesh.faces, pose, time,
        center_gravity=[0, 0, -3.9], rho=1000, gravity=gravity,
        water_depth=depth, wave_number=wave_number, wave_height=2.5,
        wave_period=8, ramp_time=40,
    )
    hydro_gap = np.abs(pressure.hydrostatic - source["hydrostatic"])
    nonlinear_gap = np.abs(pressure.nonlinear_wave - source["nonlinearWave"])
    linear_gap = np.abs(pressure.linear_wave - source["linearWave"])
    # A sparse facet can switch wet/dry classification when its center lies
    # on the moving surface. The two jump terms cancel in physical pressure.
    switches = (hydro_gap > 10) | (nonlinear_gap > 10)
    assert np.count_nonzero(switches) <= 25
    assert np.array_equal(hydro_gap > 10, nonlinear_gap > 10)
    assert np.max(hydro_gap[~switches]) < 8
    assert np.max(nonlinear_gap[~switches]) < 2
    _maximum(pressure.hydrostatic + pressure.nonlinear_wave,
             source["hydrostatic"] + source["nonlinearWave"],
             12, "combined wet-surface pressure")
    assert np.max(linear_gap) < 2

    # This is the complete 1,201-frame scene driven by the Python states and
    # pressures above, with the project's unchanged PVD/metadata writers.
    wave = WaveClass("regular")
    wave.A, wave.k, wave.w = 1.25, wave_number, omega
    wave.waveDir = [0]
    wave.waterDepth = depth
    writer = ParaviewClass(wave)
    wave_paths = writer.write_paraview_vtp_wave(
        time, tmp_path, domain_size=200, num_points_x=50, num_points_y=50,
    )
    flap_paths = writer.write_paraview_vtp(
        time, tmp_path, body_name="flap", vertices=mesh.vertices,
        faces=mesh.faces, poses=pose, cell_areas=mesh.area_faces,
        hydrostatic_pressure=pressure.hydrostatic,
        nonlinear_wave_pressure=pressure.nonlinear_wave,
        linear_wave_pressure=pressure.linear_wave,
    )
    base_mesh = trimesh.load_mesh(BASE, process=True)
    base_paths = writer.write_paraview_vtp(
        time, tmp_path, body_name="base", body_index=2,
        vertices=base_mesh.vertices, faces=base_mesh.faces,
        poses=result.bodies["base"].position,
    )
    collection = writer.write_paraview_response(
        time, tmp_path, model_name="OSWEC",
        bodies=[("flap", (1, 1, 0), 1), ("base", (1, 1, 0), 1)],
    )
    assert len(wave_paths) == len(flap_paths) == len(base_paths) == 1201
    entries = ElementTree.parse(collection).findall("./Collection/DataSet")
    assert len(entries) == 3603
    assert all((tmp_path / item.get("file")).is_file() for item in entries)

    with ZipFile(FIXTURE / "vtp_samples.zip") as saved:
        for name in ("ground.txt", "bodies.txt"):
            assert (tmp_path / name).read_bytes() == saved.read(name)
        source_collection = ElementTree.fromstring(saved.read("OSWEC.pvd"))
        assert [item.attrib for item in entries] == [
            item.attrib for item in source_collection.findall("./Collection/DataSet")
        ]
        for frame in (1, 601, 1201):
            _compare_scene_frame(frame, wave_paths[frame - 1],
                                 flap_paths[frame - 1], base_paths[frame - 1],
                                 saved.read)

    # The hosted MATLAB job retains all freshly generated source VTP frames.
    # Compare every file there; ordinary CI uses the portable numerical and
    # scene gates above plus the three pinned source VTP samples.
    source_vtp = os.environ.get("WEC_SIM_OSWEC_CORRECTED_VTP_DIR")
    if source_vtp:
        source_vtp = Path(source_vtp)
        for frame in range(1, 1202):
            _compare_scene_frame(
                frame, wave_paths[frame - 1], flap_paths[frame - 1],
                base_paths[frame - 1],
                lambda name: (source_vtp / name).read_bytes(),
            )


@pytest.mark.skipif(
    not os.environ.get("WEC_SIM_OSWEC_CORRECTED_DIR"),
    reason="fresh corrected MATLAB output not provided",
)
def test_corrected_archive_regenerates_from_pinned_matlab():
    from scipy.io import loadmat

    directory = Path(os.environ["WEC_SIM_OSWEC_CORRECTED_DIR"])
    with np.load(FIXTURE / "history.npz", allow_pickle=False) as archive:
        for filename, prefix in (("baseline.mat", ""),
                                 ("stability.mat", "stable_")):
            fresh = loadmat(directory / filename)
            for field in ("pose", "velocity", "basePose", "baseVelocity",
                          "ptoPosition", "ptoVelocity", "ptoForce", "ptoPower"):
                _maximum(fresh[field], archive[f"{prefix}{field}"],
                         1e-9, f"regenerated {prefix}{field}")
            if not prefix:
                for field in ("hydrostatic", "nonlinearWave", "linearWave"):
                    _maximum(fresh[field], archive[field], 1e-6,
                             f"regenerated {field}")

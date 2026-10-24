"""Portable 80 s RM3 MoorDyn corrected-baseline and own-state scene gate."""

import hashlib
import os
from pathlib import Path
import shutil
from xml.etree import ElementTree

import h5py
import numpy as np
import pytest
import trimesh

from wecsim import JONSWAPWave, MoorDyn, WEC
from wecsim.environment.irregularWave import jonswap_equal_energy_components
from wecsim.environment.waveClass import WaveClass
from wecsim.visualization.paraviewClass import ParaviewClass


ROOT = Path(__file__).parent.parent
FIXTURE = ROOT / "tests/fixtures/matlab/rm3_moordyn_viz_corrected"
PUBLISHED = ROOT / "tests/fixtures/matlab/rm3_moordyn_source_step.npz"
LIBRARY = os.getenv("WEC_SIM_MOORDYN_LIBRARY")
FRESH_COMPACT = os.getenv("WEC_SIM_MATLAB_RM3_VIZ_COMPACT")
APPLICATIONS = os.getenv("WEC_SIM_APPLICATIONS_DIR")

# The published gate in test_rm3_moordyn_viz_source.py retains these same
# numerical limits and its strict expected failure on the original 0.01 s
# MATLAB trajectory. This adopted baseline changes only the source coupling
# step, fixes the random seed, and keeps the original 0.1 s output/scene grid.
BODY_POSITION_LIMITS = (0.01, 0.005, 0.0005)
BODY_VELOCITY_LIMITS = (0.005, 0.005, 0.0005)
CONNECTION_LIMITS = {
    "position": (0.01, 1e-6, 0.005, 1e-6, 0.0005, 1e-6),
    "velocity": (0.005, 1e-6, 0.005, 1e-6, 0.0005, 1e-6),
    "force": (2000, 1e-3, 4000, 1e-3, 2000, 1e-3),
}
FINE_STEP_LIMITS = {"position": (0.002, 5e-5),
                    "velocity": (0.0005, 5e-5),
                    "force": (500, 2000)}
STL_SHA256 = {
    "float.stl": (
        "3d990ab8584fa21af2f22278bfe862fe119a491e4caa0d8d9abd988d1096bf4b"
    ),
    "plate.stl": (
        "9db8cc11ecf74373c0b570651e91e822617379b58f4137432028cfbb1d107acb"
    ),
}


def _maximum(actual, expected):
    assert actual.shape == expected.shape
    assert np.isfinite(actual).all() and np.isfinite(expected).all()
    return float(np.max(np.abs(actual - expected)))


def _bound(actual, expected, limit, label):
    error = _maximum(actual, expected)
    print(f"{label}: {error:.8g} (gate {limit:.8g})")
    assert np.isfinite(error) and error < limit, label


def _vtp_points(root, selector):
    value = root.findtext(selector)
    assert value is not None
    values = np.fromstring(value, sep=" ")
    assert values.size % 3 == 0 and np.isfinite(values).all()
    return values.reshape(-1, 3)


def _body_stats(path, expected_facets):
    root = ElementTree.parse(path).getroot()
    piece = root.find("./PolyData/Piece")
    points = _vtp_points(root, "./PolyData/Piece/Points/DataArray")
    connectivity = np.fromstring(piece.findtext(
        "./Polys/DataArray[@Name='connectivity']"
    ), sep=" ", dtype=int)
    offsets = np.fromstring(piece.findtext(
        "./Polys/DataArray[@Name='offsets']"
    ), sep=" ", dtype=int)
    cell_area = np.fromstring(piece.findtext(
        "./CellData/DataArray[@Name='CellArea']"
    ), sep=" ")
    facets = int(piece.get("NumberOfPolys"))
    assert (facets == expected_facets
            and len(points) == int(piece.get("NumberOfPoints"))
            and connectivity.shape == (facets * 3,)
            and offsets.shape == cell_area.shape == (facets,)
            and np.array_equal(offsets, 3 * np.arange(1, facets + 1))
            and np.all((connectivity >= 0) & (connectivity < len(points)))
            and np.isfinite(cell_area).all() and np.all(cell_area > 0))
    triangles = points[connectivity.reshape(-1, 3)]
    centers = triangles.mean(axis=1)
    area = np.linalg.norm(np.cross(
        triangles[:, 1] - triangles[:, 0],
        triangles[:, 2] - triangles[:, 0],
    ), axis=1) / 2
    return np.r_[centers.mean(axis=0), points.min(axis=0),
                 points.max(axis=0), area.sum(), cell_area.sum()]


def test_corrected_source_retains_sea_and_resolves_step_sensitivity():
    assert hashlib.sha256((FIXTURE / "rm3.h5").read_bytes()).hexdigest() == (
        "e0c163a8f2e2254e63fb6261fac8b3e539592efb8f5d4db7ec7e0f5851668070"
    )
    assert hashlib.sha256((FIXTURE / "lines.txt").read_bytes()).hexdigest() == (
        "d98ac9edcd8077363fe899f644a88ffc2f12ab88215ac2c35fdf6b6e0af99ded"
    )
    with np.load(FIXTURE / "source.npz", allow_pickle=False) as source, \
            np.load(PUBLISHED, allow_pickle=False) as published:
        assert source["phase"].shape == (3, 1000)
        np.testing.assert_array_equal(source["phase"][0], source["phase"][1])
        np.testing.assert_array_equal(source["phase"][0], source["phase"][2])
        np.testing.assert_array_equal(source["phase"][2], published["phase"][0])
        np.testing.assert_allclose(source["time"], np.arange(801) * .1,
                                   rtol=0, atol=2e-14)
        for arm in (0, 1, 2):
            _bound(source["wave"][arm], published["wave"][0], 2e-14,
                   f"arm {arm} published sea")
        for left, right, name in ((0, 1, "fine-step stability"),
                                  (1, 2, "published-output correction")):
            for kind in ("position", "velocity"):
                field = source[f"body_{kind}"]
                for body in (0, 1):
                    for axis in (0, 2, 4):
                        _bound(field[left, body, :, axis],
                               field[right, body, :, axis],
                               FINE_STEP_LIMITS[kind][axis == 4],
                               f"{name} body {body + 1} {kind} {axis}")
            for kind in ("position", "velocity", "force"):
                field = source[f"connection_{kind}"]
                for axis in (0, 2, 4):
                    _bound(field[left, :, axis], field[right, :, axis],
                           FINE_STEP_LIMITS[kind][axis == 4],
                           f"{name} connection {kind} {axis}")
            _bound(source["pto_force"][left], source["pto_force"][right],
                   2000, f"{name} PTO force")
            _bound(source["fairlead_tension"][left],
                   source["fairlead_tension"][right], 1500,
                   f"{name} fairlead tension")
            _bound(source["fairlead_dense_tension"][left],
                   source["fairlead_dense_tension"][right], 1500,
                   f"{name} 0.01 s fairlead tension")
            _bound(source["line_nodes"][left], source["line_nodes"][right],
                   .005, f"{name} line nodes")
            _bound(source["line_tension"][left],
                   source["line_tension"][right], 1500,
                   f"{name} line segment tension")
        # Keep the unchanged, strictly expected-failing published source
        # visibly separate from the corrected baseline.
        assert _maximum(source["body_position"][2, 0, :, 0],
                        published["body1_position"][0, :, 0]) > .18
        assert _maximum(source["connection_force"][2, :, 4],
                        published["mooring_force"][0, :, 4]) > 100_000
        _bound(source["scene_line_nodes"], source["line_nodes"][2],
               1.01e-5, "corrected source line VTP positions")
        _bound(source["scene_line_tension"], source["line_tension"][2],
               .1, "corrected source line VTP tensions")
        np.testing.assert_allclose(source["fairlead_dense_time"],
                                   np.arange(8001) * .01, rtol=0,
                                   atol=2e-14)
        np.testing.assert_allclose(source["fairlead_dense_tension"][:, ::10],
                                   source["fairlead_tension"], rtol=0,
                                   atol=1e-8)


@pytest.mark.skipif(not (FRESH_COMPACT and APPLICATIONS),
                    reason="fresh corrected MATLAB archive and inputs absent")
def test_fresh_corrected_source_reproduces_portable_fixture():
    apps = Path(APPLICATIONS)
    published_input = (apps / "Paraview_Visualization/RM3_MoorDyn_Viz"
                       / "wecSimInputFile.m")
    corrected_input = (apps / "Paraview_Visualization/RM3_MoorDyn_Viz_corrected"
                       / "wecSimInputFile.m")
    assert hashlib.sha256(published_input.read_bytes()).hexdigest() == (
        "16a0d1519e218429924bcd879335abdb86b51bbf91fa0f5a0a781d93a6ad9ff5"
    )
    assert hashlib.sha256(corrected_input.read_bytes()).hexdigest() == (
        "f403e48a842fb1c0fe73b73d7476dcf4b652794a449bee08a077f3c66dfe8df1"
    )
    source_lines = (apps / "Paraview_Visualization/RM3_MoorDyn_Viz"
                    / "Mooring/lines.txt")
    assert source_lines.read_bytes() == (FIXTURE / "lines.txt").read_bytes()
    live_hydro = apps / "_Common_Input_Files/RM3/hydroData/rm3.h5"
    with h5py.File(live_hydro) as live, h5py.File(FIXTURE / "rm3.h5") as saved:
        def datasets(group):
            found = {}
            group.visititems(lambda name, item: found.__setitem__(
                name, item[()]
            ) if isinstance(item, h5py.Dataset) else None)
            return found

        live_fields, saved_fields = datasets(live), datasets(saved)
        assert live_fields.keys() == saved_fields.keys()
        for name, values in live_fields.items():
            reference = saved_fields[name]
            assert values.shape == reference.shape, name
            if np.issubdtype(values.dtype, np.number):
                if ("/excitation/froude-krylov/" not in name and
                        "/excitation/scattering/" not in name):
                    assert (np.isfinite(values).all() and
                            np.isfinite(reference).all()), name
            # Regenerated HDF5 containers can differ bytewise; all 84
            # datasets, including active FFT IRF K, must be exact. BEMIO's
            # unused Froude-Krylov/scattering fields contain matching NaNs.
            np.testing.assert_array_equal(values, reference,
                                          err_msg=f"fresh HDF5 {name}")
    for stl, digest in STL_SHA256.items():
        source_mesh = apps / "_Common_Input_Files/RM3/geometry" / stl
        assert hashlib.sha256(source_mesh.read_bytes()).hexdigest() == digest
        assert hashlib.sha256((ROOT / "examples/data" / stl).read_bytes()).hexdigest() == digest

    limits = {
        "time": 1e-8,
        "phase": 1e-12,
        "wave": 1e-10,
        "body_position": 1e-7,
        "body_velocity": 1e-7,
        "pto_force": .1,
        "pto_power": .1,
        "connection_position": 1e-7,
        "connection_velocity": 1e-7,
        "connection_force": .1,
        "fairlead_tension": .1,
        "fairlead_dense_time": 1e-8,
        "fairlead_dense_tension": .1,
        "line_nodes": 1e-7,
        "line_tension": .1,
        "scene_wave_samples": 1.01e-5,
        "scene_wave_z_extrema": 1.01e-5,
        "scene_line_nodes": 1.01e-5,
        "scene_line_tension": .1,
    }
    with np.load(FIXTURE / "source.npz", allow_pickle=False) as saved, \
            np.load(FRESH_COMPACT, allow_pickle=False) as fresh:
        assert set(saved.files) == set(fresh.files) == (
            set(limits) | {"scene_wave_connectivity",
                           "scene_wave_offsets", "scene_wave_sample_indices",
                           "scene_body_stats"}
        )
        for key in saved.files:
            if key == "scene_body_stats":
                _bound(fresh[key][..., :9], saved[key][..., :9],
                       1.01e-5, "fresh archive body scene coordinates")
                _bound(fresh[key][..., 9:], saved[key][..., 9:],
                       1e-3, "fresh archive body scene areas")
            elif key in limits:
                _bound(fresh[key], saved[key], limits[key],
                       f"fresh archive {key}")
            else:
                np.testing.assert_array_equal(fresh[key], saved[key])


@pytest.mark.skipif(not LIBRARY, reason="pinned native MoorDyn library absent")
def test_independent_full_duration_motion_and_own_state_scene(tmp_path):
    for stl, digest in STL_SHA256.items():
        assert (hashlib.sha256((ROOT / "examples/data" / stl).read_bytes())
                .hexdigest() == digest)
    with np.load(FIXTURE / "source.npz", allow_pickle=False) as source:
        expected = {name: source[name].copy() for name in source.files}
    phase_file = tmp_path / "phase.csv"
    np.savetxt(phase_file, expected["phase"][2], delimiter=",")
    mooring_input = tmp_path / "Mooring" / "lines.txt"
    mooring_input.parent.mkdir()
    shutil.copyfile(FIXTURE / "lines.txt", mooring_input)
    wec = WEC("Corrected RM3 MoorDyn visualization")
    hydro = FIXTURE / "rm3.h5"
    float_body = wec.body("float", hydro, inertia=(0, 21_306_090.66, 0))
    spar = wec.body("spar", hydro, inertia=(0, 94_407_091.24, 0))
    moordyn = MoorDyn(LIBRARY, mooring_input)
    wec.floating_joint(
        float_body, spar, damping=1_200_000, moordyn=moordyn,
        moordyn_point=spar.at(0, 0, 21.5),
    )
    result = wec.run(
        JONSWAPWave(2, 8, phase_file=phase_file,
                    discretization="traditional"),
        dt=.00125, end_time=80, ramp_time=0, radiation_memory=60,
        initial_coordinate={"spar_heave": -.21},
    )
    assert result.time.shape == (64001,)
    indices = np.arange(801) * 80
    _bound(result.time[indices], expected["time"], 1e-8,
           "independent sample time")
    _bound(result.wave_elevation[indices], expected["wave"][2],
           1e-10, "independent incident sea")
    for body, name in enumerate(("float", "spar")):
        actual = result.bodies[name]
        for kind, limits in (("position", BODY_POSITION_LIMITS),
                             ("velocity", BODY_VELOCITY_LIMITS)):
            for axis, limit in zip((0, 2, 4), limits, strict=True):
                _bound(getattr(actual, kind)[indices, axis],
                       expected[f"body_{kind}"][2, body, :, axis], limit,
                       f"{name} {kind} {axis}")
    pto = result.ptos["relative_heave"]
    _bound(pto.force[indices], expected["pto_force"][2], 2000,
           "independent PTO force")
    for start, stop in ((0, 80), (0, 20), (20, 50), (50, 80)):
        window = (expected["time"] >= start) & (expected["time"] <= stop)
        source_work = np.trapezoid(-expected["pto_power"][2, window],
                                   expected["time"][window])
        python_work = np.trapezoid(pto.absorbed_power[indices][window],
                                   expected["time"][window])
        assert source_work > 1e6
        work_error = abs(python_work - source_work) / source_work
        print(f"PTO signed work {start:g}-{stop:g} s: {work_error:.8g} "
              "(gate 0.01)")
        assert work_error < .01
    outputs = dict(result.raw.extra_outputs)
    for kind, limits in CONNECTION_LIMITS.items():
        actual = outputs[f"moordyn_connection_{kind}"][indices]
        for axis, limit in enumerate(limits):
            _bound(actual[:, axis],
                   expected[f"connection_{kind}"][2, :, axis], limit,
                   f"independent connection {kind} {axis}")
    fairlead = np.loadtxt(mooring_input.with_suffix(".out"),
                          skiprows=1, ndmin=2)
    assert fairlead.shape == (64001, 4)
    assert np.isfinite(fairlead).all()
    dense_indices = np.arange(8001) * 8
    _bound(fairlead[dense_indices, 0],
           expected["fairlead_dense_time"], 1e-8,
           "native 0.01 s fairlead time")
    _bound(fairlead[dense_indices, 1:4],
           expected["fairlead_dense_tension"][2], 4000,
           "independent 0.01 s fairlead tensions")
    _bound(fairlead[indices, 0], expected["time"], 1e-8,
           "native fairlead time")
    _bound(fairlead[indices, 1:4], expected["fairlead_tension"][2],
           4000, "independent fairlead tensions")

    native_nodes, native_tensions = [], []
    for line in (1, 2, 3):
        time, nodes, tension = moordyn.line_history(line)
        assert nodes.shape == (64001, 21, 3)
        assert tension.shape == (64001, 20)
        assert np.isfinite(time).all()
        assert np.isfinite(nodes).all() and np.isfinite(tension).all()
        _bound(time[indices], expected["time"], 1e-8,
               f"native line {line} frame time")
        _bound(nodes[indices], expected["line_nodes"][2, line - 1],
               .01, f"independent line {line} nodes")
        _bound(tension[indices], expected["line_tension"][2, line - 1],
               4000, f"independent line {line} segment tensions")
        native_nodes.append(nodes[indices])
        native_tensions.append(tension[indices])

    # Write every published 0.1 s visualization frame from this independent
    # trajectory, then compare the actual XML output to the corrected MATLAB
    # scene. No saved MATLAB body, connection, or line state drives a writer.
    scene_dir = tmp_path / "own_scene"
    paths = ParaviewClass(None).write_paraview_vtp_mooring(
        expected["time"], expected["time"], scene_dir,
        lines=native_nodes, tensions=native_tensions,
    )
    assert len(paths) == 801
    scene_node_error = 0.0
    scene_tension_error = 0.0
    for frame, path in enumerate(paths):
        pieces = ElementTree.parse(path).findall("./PolyData/Piece")
        assert len(pieces) == 3
        for line, piece in enumerate(pieces):
            assert piece.get("NumberOfPoints") == "21"
            assert piece.get("NumberOfLines") == "20"
            points = _vtp_points(piece, "./Points/DataArray")
            connectivity = np.fromstring(piece.findtext(
                "./Lines/DataArray[@Name='connectivity']"
            ), sep=" ", dtype=int).reshape(-1, 2)
            offsets = np.fromstring(piece.findtext(
                "./Lines/DataArray[@Name='offsets']"
            ), sep=" ", dtype=int)
            load = np.fromstring(piece.findtext(
                "./CellData/DataArray[@Name='Segment Tension']"
            ), sep=" ")
            assert load.shape == (20,) and np.isfinite(load).all()
            np.testing.assert_array_equal(connectivity, np.column_stack((
                np.arange(20), np.arange(1, 21),
            )))
            np.testing.assert_array_equal(offsets, 2 * np.arange(1, 21))
            scene_node_error = max(scene_node_error, _maximum(
                points, expected["scene_line_nodes"][line, frame],
            ))
            scene_tension_error = max(scene_tension_error, _maximum(
                load, expected["scene_line_tension"][line, frame],
            ))
    assert scene_node_error < .01002
    assert scene_tension_error < 4000.1
    print(f"own-scene line VTP nodes: {scene_node_error:.8g} "
          "(gate 0.01002)")
    print(f"own-scene line VTP tensions: {scene_tension_error:.8g} "
          "(gate 4000.1)")

    for body, (name, stl) in enumerate((("float", "float.stl"),
                                         ("spar", "plate.stl")), start=1):
        mesh = trimesh.load_mesh(ROOT / "examples/data" / stl,
                                 process=False)
        paths = ParaviewClass(None).write_paraview_vtp(
            expected["time"], scene_dir, body_name=name,
            vertices=mesh.vertices, faces=mesh.faces,
            poses=result.bodies[name].position[indices], body_index=body,
        )
        assert len(paths) == 801
        assert len(mesh.faces) == (516, 356)[body - 1]
        stats = np.stack([_body_stats(path, len(mesh.faces)) for path in paths])
        _bound(stats[:, :9], expected["scene_body_stats"][body - 1, :, :9],
               .005, f"own-scene {name} geometry")
        _bound(stats[:, 9], expected["scene_body_stats"][body - 1, :, 9],
               .005, f"own-scene {name} area")
        _bound(stats[:, 10], expected["scene_body_stats"][body - 1, :, 10],
               .005, f"own-scene {name} CellArea")

    components = jonswap_equal_energy_components(
        hydro, significant_height=2, peak_period=8,
        directions=[0], spreading=[1], phase=expected["phase"][2, :, None],
        discretization="traditional",
    )
    wave = WaveClass("irregular")
    wave.w = components.omega
    # The source VTP wave.A field is twice its archived component density;
    # jonswap_equal_energy_components already returns that VTP convention.
    wave.A = components.spectral_amplitude
    wave.dw = components.d_omega
    wave.phase = components.phase.T
    wave.waveDir = [0]
    wave.waveSpread = [1]
    wave.waterDepth = 70
    wave.deepWaterWave = 0
    wave.waveNumber(9.81)
    wave.viz = {"numPointsX": 1000, "numPointsY": 2}
    paths = ParaviewClass(wave).write_paraview_vtp_wave(
        expected["time"], scene_dir, domain_size=300, num_moordyn=1,
    )
    assert len(paths) == 801
    assert (scene_dir / "ground.txt").read_text() == "300\n70\n1\n"
    scene_wave_sample_error = 0.0
    scene_wave_extent_error = 0.0
    for frame, path in enumerate(paths):
        root = ElementTree.parse(path).getroot()
        points = _vtp_points(root,
                             "./PolyData/Piece/Points/DataArray")
        assert points.shape == (2000, 3)
        piece = root.find("./PolyData/Piece")
        assert piece.get("NumberOfPoints") == "2000"
        assert piece.get("NumberOfPolys") == "999"
        connectivity = np.fromstring(piece.findtext(
            "./Polys/DataArray[@Name='connectivity']"
        ), sep=" ", dtype=int)
        offsets = np.fromstring(piece.findtext(
            "./Polys/DataArray[@Name='offsets']"
        ), sep=" ", dtype=int)
        np.testing.assert_array_equal(connectivity,
                                      expected["scene_wave_connectivity"])
        np.testing.assert_array_equal(offsets, 4 * np.arange(1, 1000))
        scene_wave_sample_error = max(scene_wave_sample_error, _maximum(
            points[expected["scene_wave_sample_indices"]],
            expected["scene_wave_samples"][frame],
        ))
        scene_wave_extent_error = max(scene_wave_extent_error, _maximum(
            np.array([points[:, 2].min(), points[:, 2].max()]),
            expected["scene_wave_z_extrema"][frame],
        ))
    assert scene_wave_sample_error < 1.01e-5
    assert scene_wave_extent_error < 1.01e-5
    print(f"own-scene wave sampled grid: {scene_wave_sample_error:.8g} "
          "(gate 1.01e-5)")
    print(f"own-scene wave z extrema: {scene_wave_extent_error:.8g} "
          "(gate 1.01e-5)")

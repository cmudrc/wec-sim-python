"""Pair the published Traditional JONSWAP RM3 MoorDyn visualization motion."""

import os
from pathlib import Path
import shutil
from xml.etree import ElementTree

import h5py
import numpy as np
import pytest
import trimesh
from scipy.spatial import cKDTree

from wecsim import JONSWAPWave, MoorDyn, WEC
from wecsim.irregularWave import (
    jonswap_equal_energy_components, synthesize_irregular_response,
)
from wecsim.paraviewClass import ParaviewClass
from wecsim.waveClass import WaveClass


APPLICATIONS = os.environ.get("WEC_SIM_APPLICATIONS_DIR")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_RM3_VIZ_DIR")
LIBRARY = os.environ.get("WEC_SIM_MOORDYN_LIBRARY")
pytestmark = pytest.mark.skipif(
    not (APPLICATIONS and REFERENCE and LIBRARY),
    reason="pinned RM3 MoorDyn visualization input, output, and library not provided",
)


def _read(name):
    values = np.loadtxt(Path(REFERENCE) / name, delimiter=",", ndmin=2)
    assert np.isfinite(values).all(), name
    return values


def _bound(actual, expected, limit, name):
    assert actual.shape == expected.shape, name
    error = float(np.max(np.abs(actual - expected)))
    print(f"{name}: maximum difference {error:.8g}; gate {limit:.8g}")
    assert error < limit, name
    return error


FINE_SOURCE_LIMITS = {
    "position": (0.002, 5e-5),
    "velocity": (0.0005, 5e-5),
    "force": (500, 2000),
}


class PublishedCoupledMotionGap(AssertionError):
    """The published coarse-step trajectory misses its unchanged paired gates."""


EXPECTED_PUBLISHED_GAPS = {
    "float position axis 0", "float velocity axis 0",
    "float position axis 4", "float velocity axis 4",
    "spar position axis 0", "spar velocity axis 0",
    "spar position axis 4", "spar velocity axis 4",
    "mooring position axis 0", "mooring position axis 4",
    "mooring velocity axis 0", "mooring velocity axis 4",
    "mooring force axis 0", "mooring force axis 2",
    "mooring force axis 4",
    "three fairlead tensions",
}


def _indices(source_time, *, dt=0.01, end_time=80):
    indices = np.rint(source_time / dt).astype(int)
    assert np.all((indices >= 0) & (indices <= round(end_time / dt)))
    np.testing.assert_allclose(indices * dt, source_time, rtol=0, atol=1e-8)
    return indices


def _wave_vtp(path):
    piece = ElementTree.parse(path).find("./PolyData/Piece")
    points = np.fromstring(piece.findtext("./Points/DataArray"), sep=" ")
    cells = np.fromstring(
        piece.findtext("./Polys/DataArray[@Name='connectivity']"),
        sep=" ", dtype=int,
    )
    offsets = np.fromstring(
        piece.findtext("./Polys/DataArray[@Name='offsets']"),
        sep=" ", dtype=int,
    )
    return points.reshape(-1, 3), cells.reshape(-1, 4), offsets


def _body_vtp(path):
    piece = ElementTree.parse(path).find("./PolyData/Piece")
    points = np.fromstring(piece.findtext("./Points/DataArray"), sep=" ").reshape(-1, 3)
    faces = np.fromstring(
        piece.findtext("./Polys/DataArray[@Name='connectivity']"),
        sep=" ", dtype=int,
    ).reshape(-1, 3)
    assert len(points) == int(piece.get("NumberOfPoints"))
    assert len(faces) == int(piece.get("NumberOfPolys"))
    cell_data = piece.findall("./CellData/DataArray")
    assert [field.get("Name") for field in cell_data] == ["CellArea"]
    cell_area = np.fromstring(cell_data[0].text, sep=" ")
    assert cell_area.shape == (len(faces),)
    return points[faces], cell_area


def _mooring_vtp(path):
    pieces = ElementTree.parse(path).findall("./PolyData/Piece")
    assert len(pieces) == 3
    records = []
    for piece in pieces:
        points = np.fromstring(
            piece.findtext("./Points/DataArray"), sep=" ",
        ).reshape(-1, 3)
        cells = np.fromstring(
            piece.findtext("./Lines/DataArray[@Name='connectivity']"),
            sep=" ", dtype=int,
        ).reshape(-1, 2)
        offsets = np.fromstring(
            piece.findtext("./Lines/DataArray[@Name='offsets']"),
            sep=" ", dtype=int,
        )
        tension = np.fromstring(
            piece.findtext("./CellData/DataArray[@Name='Segment Tension']"),
            sep=" ",
        )
        assert (points.shape == (21, 3) and cells.shape == (20, 2)
                and offsets.shape == tension.shape == (20,))
        records.append((points, cells, offsets, tension))
    return records


def test_published_moordyn_line_vtp_history(tmp_path):
    """Pair every actual line frame from saved source nodes and tensions."""
    source = Path(REFERENCE)
    histories = [_read(f"line{line}_vtp_history.csv") for line in (1, 2, 3)]
    assert all(history.shape == (801, 84) for history in histories)
    times = histories[0][:, 0]
    np.testing.assert_allclose(times, np.arange(801) * 0.1, rtol=0, atol=1e-8)
    for history in histories[1:]:
        np.testing.assert_allclose(history[:, 0], times, rtol=0, atol=1e-8)
    lines = [history[:, 1:64].reshape(801, 21, 3) for history in histories]
    tensions = [history[:, 64:] for history in histories]
    paths = ParaviewClass(None).write_paraview_vtp_mooring(
        times, times, tmp_path, lines=lines, tensions=tensions,
    )
    assert len(paths) == 801
    max_point_error = 0.0
    max_tension_error = 0.0
    for frame, path in enumerate(paths, start=1):
        expected = _mooring_vtp(
            source / "published_vtp/mooring1" / f"mooring_{frame}.vtp"
        )
        actual = _mooring_vtp(path)
        for actual_line, expected_line in zip(actual, expected):
            points, cells, offsets, load = actual_line
            source_points, source_cells, source_offsets, source_load = expected_line
            max_point_error = max(
                max_point_error, float(np.max(np.abs(points - source_points))),
            )
            max_tension_error = max(
                max_tension_error, float(np.max(np.abs(load - source_load))),
            )
            np.testing.assert_array_equal(cells, source_cells)
            np.testing.assert_array_equal(offsets, source_offsets)
    print(f"published MoorDyn VTP node difference: {max_point_error:.8g} m")
    print(f"published MoorDyn VTP tension difference: {max_tension_error:.8g} N")
    assert max_point_error < 1.01e-5
    assert max_tension_error < 0.1


def test_published_body_vtp_history_on_source_poses(tmp_path):
    """Pair every body mesh frame; source poses are writer inputs."""
    reference = Path(REFERENCE)
    apps = Path(APPLICATIONS)
    for body_index, (name, stl) in enumerate(
        (("float", "float.stl"), ("spar", "plate.stl")), start=1,
    ):
        states = _read(f"body{body_index}.csv")
        assert states.shape == (801, 25)
        np.testing.assert_allclose(states[:, 0], np.arange(801) * .1,
                                   rtol=0, atol=1e-8)
        mesh = trimesh.load_mesh(
            apps / "_Common_Input_Files/RM3/geometry" / stl, process=False,
        )
        paths = ParaviewClass(None).write_paraview_vtp(
            states[:, 0], tmp_path, body_name=name,
            vertices=mesh.vertices, faces=mesh.faces,
            poses=states[:, 1:7], body_index=body_index,
        )
        assert len(paths) == 801
        max_center_error = 0.0
        max_vertex_error = 0.0
        max_area_error = 0.0
        max_cell_area_error = 0.0
        for frame, path in enumerate(paths, start=1):
            source = (reference / "published_vtp" / f"body{body_index}_{name}"
                      / f"{name}_{frame}.vtp")
            expected, source_cell_area = _body_vtp(source)
            actual, actual_cell_area = _body_vtp(path)
            assert actual.shape == expected.shape
            # STL readers can number shared vertices differently. Compare
            # every geometric facet using a one-to-one centroid match.
            source_center = expected.mean(axis=1)
            actual_center = actual.mean(axis=1)
            distance, match = cKDTree(actual_center).query(source_center)
            assert len(np.unique(match)) == len(match)
            max_center_error = max(max_center_error, float(np.max(distance)))
            # Compare all three vertices of every matched facet, independent
            # of each STL reader's winding and starting-vertex convention.
            vertex_distance = np.linalg.norm(
                expected[:, :, None, :] - actual[match][:, None, :, :],
                axis=3,
            )
            max_vertex_error = max(
                max_vertex_error,
                float(np.max(np.min(vertex_distance, axis=2))),
                float(np.max(np.min(vertex_distance, axis=1))),
            )
            source_area = np.linalg.norm(np.cross(
                expected[:, 1] - expected[:, 0],
                expected[:, 2] - expected[:, 0],
            ), axis=1) / 2
            actual_area = np.linalg.norm(np.cross(
                actual[:, 1] - actual[:, 0],
                actual[:, 2] - actual[:, 0],
            ), axis=1) / 2
            max_area_error = max(max_area_error, float(np.max(
                np.abs(actual_area[match] - source_area),
            )))
            max_cell_area_error = max(max_cell_area_error, float(np.max(
                np.abs(actual_cell_area[match] - source_cell_area),
            )))
        print(f"published {name} VTP maximum facet-center difference: "
              f"{max_center_error:.8g} m")
        print(f"published {name} VTP maximum facet-vertex difference: "
              f"{max_vertex_error:.8g} m")
        print(f"published {name} VTP maximum facet-area difference: "
              f"{max_area_error:.8g} m²")
        print(f"published {name} VTP maximum CellArea difference: "
              f"{max_cell_area_error:.8g} m²")
        assert max_center_error < 3e-5
        assert max_vertex_error < 3e-5
        assert max_area_error < 1e-3
        assert max_cell_area_error < 2e-5


def test_published_irregular_wave_vtp_history(tmp_path):
    """Pair every actual RM3 wave mesh from the 80 s application."""
    reference = Path(REFERENCE)
    parameters = _read("wave_vtp_parameters.csv").ravel()
    assert parameters.shape == (7,)
    domain, depth, moorings, nx, ny, direction, spreading = parameters
    assert (domain, depth, moorings, nx, ny) == (300, 70, 1, 1000, 2)
    components = _read("wave_vtp_components.csv")
    assert components.shape == (1000, 5)
    sea_components = _read("components.csv")
    _bound(components[:, (0, 3, 4)], sea_components[:, (0, 1, 3)],
           1e-12, "published VTP frequency, width, and phase")
    _bound(components[:, 2], 2 * sea_components[:, 2],
           1e-12, "published VTP amplitude convention")
    wave = WaveClass("irregular")
    wave.w = components[:, 0]
    wave.k = components[:, 1]
    wave.A = components[:, 2]
    wave.dw = components[:, 3]
    wave.phase = components[:, 4][None, :]
    wave.waveDir = [direction]
    wave.waveSpread = [spreading]
    wave.waterDepth = depth
    wave.viz = {"numPointsX": int(nx), "numPointsY": int(ny)}
    source_wave = _read("wave.csv")
    assert source_wave.shape == (801, 2)
    np.testing.assert_allclose(source_wave[:, 0], np.arange(801) * .1,
                               rtol=0, atol=1e-8)
    origin = np.array([
        ParaviewClass(wave).waveElevationGrid(
            time, np.array([[0.]]), np.array([[0.]]),
        ).item()
        for time in (0, 10, 80)
    ])
    _bound(origin, source_wave[[0, 100, 800], 1],
           1e-10, "published VTP component origin elevation")
    paths = ParaviewClass(wave).write_paraview_vtp_wave(
        source_wave[:, 0], tmp_path, domain_size=domain,
        num_moordyn=int(moorings),
    )
    assert len(paths) == 801
    assert (tmp_path / "ground.txt").read_text() == (
        reference / "published_vtp/ground.txt"
    ).read_text()
    max_point_error = 0.0
    for frame, path in enumerate(paths, start=1):
        actual_points, actual_cells, actual_offsets = _wave_vtp(path)
        source_points, source_cells, source_offsets = _wave_vtp(
            reference / "published_vtp/waves" / f"waves_{frame}.vtp"
        )
        assert actual_points.shape == source_points.shape == (2000, 3)
        assert actual_cells.shape == source_cells.shape == (999, 4)
        max_point_error = max(max_point_error, float(np.max(
            np.abs(actual_points - source_points),
        )))
        np.testing.assert_array_equal(actual_cells, source_cells)
        np.testing.assert_array_equal(actual_offsets, source_offsets)
    print(f"published wave VTP maximum point difference: {max_point_error:.8g} m")
    assert max_point_error < 1.01e-5


def test_published_traditional_sea_and_body_excitation():
    apps = Path(APPLICATIONS)
    hydro = apps / "_Common_Input_Files/RM3/hydroData/rm3.h5"
    source = _read("components.csv")
    assert source.shape == (1000, 4)
    assert np.all(np.diff(source[:, 0]) > 0)
    assert np.max(np.abs(source[:, 2])) > 0
    components = jonswap_equal_energy_components(
        hydro, significant_height=2, peak_period=8,
        directions=[0], spreading=[1], phase=source[:, 3:4],
        discretization="traditional",
    )
    generated = np.column_stack((
        components.omega, components.d_omega,
        components.spectral_amplitude / 2, components.phase[:, 0],
    ))
    np.testing.assert_allclose(generated, source, rtol=2e-12, atol=1e-14)
    wave = _read("wave.csv")
    assert wave.shape[1] == 2
    sea = tuple(synthesize_irregular_response(
        hydro, components, dt=0.01, end_time=80, ramp_time=0,
        body_number=number,
    ) for number in (1, 2))
    _bound(sea[0].elevation[_indices(wave[:, 0])], wave[:, 1],
           1e-10, "wave elevation")
    for number, incident in enumerate(sea, start=1):
        body = _read(f"body{number}.csv")
        assert body.shape[1] == 25
        _bound(incident.excitation_force[_indices(body[:, 0])],
               body[:, 19:25], 1e-4, f"body{number} excitation")


def test_saved_source_pose_through_pinned_moordyn(tmp_path):
    """Diagnose whether source poses recreate source mooring forces directly."""
    source = _read("dense_mooring.csv")
    assert source.shape == (1001, 19)
    np.testing.assert_allclose(source[:, 0], np.arange(1001) * 0.01,
                               rtol=0, atol=1e-8)
    lines_source = (Path(APPLICATIONS) /
                    "Paraview_Visualization/RM3_MoorDyn_Viz/Mooring/lines.txt")
    for mode in ("at_step_end", "at_step_start"):
        run_dir = tmp_path / mode
        run_dir.mkdir()
        lines = run_dir / "lines.txt"
        shutil.copyfile(lines_source, lines)
        actual = np.zeros((len(source), 6))
        with MoorDyn(LIBRARY, lines).start(source[0, 1:7], source[0, 7:13]) as moordyn:
            for index in range(1, len(source)):
                pose_index = index if mode == "at_step_end" else index - 1
                actual[index] = moordyn.step(
                    source[pose_index, 1:7], source[pose_index, 7:13],
                    source[index - 1, 0], 0.01,
                )
        assert np.isfinite(actual).all()
        error = actual[1:] - source[1:, 13:19]
        for axis, label in ((0, "surge force"), (2, "heave force"),
                            (4, "pitch moment")):
            print(f"source-pose MoorDyn {mode} {label}: "
                  f"maximum {np.max(np.abs(error[:, axis])):.8g}, "
                  f"RMS {np.sqrt(np.mean(error[:, axis] ** 2)):.8g}")
        if mode == "at_step_end":
            _bound(actual[1:], source[1:, 13:19], 0.01,
                   "MoorDyn on saved source connection poses")


def test_source_output_interval_sensitivity():
    """Compare the same seeded source sea at 0.01 and 0.1 s output steps."""
    _bound(_read("dense_phase.csv"), _read("coarse_phase.csv"),
           1e-12, "source phase at both output steps")
    dense_wave = _read("dense_wave.csv")
    coarse_wave = _read("coarse_wave.csv")
    assert dense_wave.shape == (1001, 2)
    assert coarse_wave.shape == (101, 2)
    _bound(dense_wave[::10], coarse_wave, 1e-10,
           "source wave at both output steps")
    for number, name in ((1, "float"), (2, "spar")):
        dense = _read(f"dense_body{number}.csv")[::10]
        coarse = _read(f"coarse_body{number}.csv")
        assert dense.shape == coarse.shape == (101, 25)
        _bound(dense[:, 1:13], coarse[:, 1:13], 1e-9,
               f"source {name} motion at both output steps")
        for axis, label in ((0, "surge"), (2, "heave"), (4, "pitch")):
            error = dense[:, 1 + axis] - coarse[:, 1 + axis]
            print(f"source output-step {name} {label}: "
                  f"maximum {np.max(np.abs(error)):.8g}, "
                  f"at 1 s {error[10]:.8g}, at 10 s {error[-1]:.8g}")
    dense_mooring = _read("dense_mooring.csv")[::10]
    coarse_mooring = _read("coarse_mooring.csv")
    assert dense_mooring.shape == coarse_mooring.shape == (101, 19)
    _bound(dense_mooring[:, 13:19], coarse_mooring[:, 13:19], 1e-6,
           "source MoorDyn loads at both output steps")
    for axis, label in ((0, "surge force"), (2, "heave force"),
                        (4, "pitch moment")):
        error = dense_mooring[:, 13 + axis] - coarse_mooring[:, 13 + axis]
        print(f"source output-step mooring {label}: "
              f"maximum {np.max(np.abs(error)):.8g}, "
              f"at 1 s {error[10]:.8g}, at 10 s {error[-1]:.8g}")


def test_seeded_full_case_shares_short_source_prefix():
    """Keep the 80 s published settings aligned with the short source audit."""
    _bound(_read("phase.csv"), _read("dense_phase.csv"),
           1e-12, "seeded full-case source phase")
    _bound(_read("wave.csv")[:101], _read("coarse_wave.csv"),
           1e-10, "seeded full-case source wave through 10 s")
    for number, name in ((1, "float"), (2, "spar")):
        _bound(_read(f"body{number}.csv")[:101, :13],
               _read(f"coarse_body{number}.csv")[:, :13],
               1e-6, f"seeded full-case source {name} motion through 10 s")
    _bound(_read("mooring.csv")[:101, 13:19],
           _read("coarse_mooring.csv")[:, 13:19],
           1e-3, "seeded full-case source mooring through 10 s")


@pytest.mark.parametrize("prefix", ["step005", "step0025", "step00125",
                                    "step000625", "step0003125"])
def test_source_maximum_step_refinement_preserves_sea(prefix):
    """Keep one sea while refining source MaxStep and MoorDyn coupling."""
    _bound(_read(f"{prefix}_phase.csv"), _read("dense_phase.csv"),
           1e-12, f"{prefix} source phase")
    wave = _read(f"{prefix}_wave.csv")
    assert wave.shape == (1001, 2)
    _bound(wave, _read("dense_wave.csv"), 1e-10,
           f"{prefix} source wave")
    for number, name in ((1, "float"), (2, "spar")):
        refined = _read(f"{prefix}_body{number}.csv")
        dense = _read(f"dense_body{number}.csv")
        assert refined.shape == dense.shape == (1001, 25)
        _bound(refined[:, 0], dense[:, 0], 1e-9,
               f"{prefix} source {name} time")
        delta = refined[:, 1:13] - dense[:, 1:13]
        print(f"{prefix} source {name} surge at 10 s: "
              f"{delta[-1, 0]:.8g} m")
    mooring = _read(f"{prefix}_mooring.csv")
    assert mooring.shape == (1001, 19)
    _bound(mooring[:, 0], _read("dense_mooring.csv")[:, 0],
           1e-9, f"{prefix} source MoorDyn time")


def test_source_refined_motion_change_diagnostic():
    """Report source motion and mooring changes across coupling steps."""
    for number, name in ((1, "float"), (2, "spar")):
        records = [_read(f"{prefix}_body{number}.csv") for prefix in
                   ("step005", "step0025", "step00125", "step000625",
                    "step0003125")]
        assert all(record.shape == (1001, 25) for record in records)
        for column, unit, channel in ((1, "m", "surge"),
                                      (3, "m", "heave"),
                                      (5, "rad", "pitch"),
                                      (7, "m/s", "surge speed"),
                                      (9, "m/s", "heave speed"),
                                      (11, "rad/s", "pitch speed")):
            for left, right, step in zip(records[:-1], records[1:],
                                         ("0.005->0.0025", "0.0025->0.00125",
                                          "0.00125->0.000625",
                                          "0.000625->0.0003125")):
                delta = right[:, column] - left[:, column]
                print(f"source {name} {channel} {step} step-change {unit}: "
                      f"max {np.max(np.abs(delta)):.8g}, "
                      f"final {delta[-1]:.8g}")
    moorings = [_read(f"{prefix}_mooring.csv") for prefix in
                ("step005", "step0025", "step00125", "step000625",
                 "step0003125")]
    assert all(record.shape == (1001, 19) for record in moorings)
    for left, right, step in zip(moorings[:-1], moorings[1:],
                                 ("0.005->0.0025", "0.0025->0.00125",
                                  "0.00125->0.000625",
                                  "0.000625->0.0003125")):
        delta = right[:, 17] - left[:, 17]
        print(f"source MoorDyn pitch moment {step} step-change N m: "
              f"max {np.max(np.abs(delta)):.8g}, final {delta[-1]:.8g}")


def test_source_finest_coupling_steps_are_within_paired_gates():
    """Bound the remaining MATLAB step change before using the fine pair."""
    for number, name in ((1, "float"), (2, "spar")):
        coarse = _read(f"step000625_body{number}.csv")
        fine = _read(f"step0003125_body{number}.csv")
        _bound(coarse[:, 0], fine[:, 0], 1e-9,
               f"source finest {name} time")
        for kind, start in (("position", 1), ("velocity", 7)):
            for axis, label in ((0, "surge"), (2, "heave"), (4, "pitch")):
                _bound(coarse[:, start + axis], fine[:, start + axis],
                       FINE_SOURCE_LIMITS[kind][axis == 4],
                       f"source finest {name} {label} {kind}")
    coarse = _read("step000625_mooring.csv")
    fine = _read("step0003125_mooring.csv")
    _bound(coarse[:, 0], fine[:, 0], 1e-9, "source finest MoorDyn time")
    for kind, start in (("position", 1), ("velocity", 7), ("force", 13)):
        for axis, label in ((0, "surge"), (2, "heave"), (4, "pitch")):
            _bound(coarse[:, start + axis], fine[:, start + axis],
                   FINE_SOURCE_LIMITS[kind][axis == 4],
                   f"source finest MoorDyn {label} {kind}")


@pytest.mark.parametrize("prefix", ["dense", "seed2_dense",
                                    "step005", "step0025", "step00125",
                                    "step000625", "step0003125"])
def test_dense_source_wave_and_force_components(prefix):
    """Check the signed WEC-Sim logging convention before using its forces."""
    hydro = (Path(APPLICATIONS) /
             "_Common_Input_Files/RM3/hydroData/rm3.h5")
    phase = _read(f"{prefix}_phase.csv")
    assert phase.shape == (1000, 1)
    seed = 2 if prefix == "seed2_dense" else 1
    components = jonswap_equal_energy_components(
        hydro, significant_height=2, peak_period=8,
        directions=[0], spreading=[1], seed=seed,
        phase_generator="matlab",
        discretization="traditional",
    )
    _bound(components.phase, phase, 1e-12,
           f"{prefix} independently generated MATLAB phase")
    wave = _read(f"{prefix}_wave.csv")
    incident = tuple(synthesize_irregular_response(
        hydro, components, dt=0.01, end_time=10, ramp_time=0,
        body_number=number,
    ) for number in (1, 2))
    _bound(incident[0].elevation, wave[:, 1], 1e-10,
           f"{prefix} source wave")
    for number in (1, 2):
        terms = _read(f"{prefix}_forces_body{number}.csv")
        assert terms.shape == (1001, 49)
        _bound(incident[number - 1].excitation_force, terms[:, 7:13],
               1e-4, f"{prefix} source body{number} excitation")
        excitation = terms[:, 7:13]
        resisting = sum(terms[:, start:start + 6]
                        for start in (13, 19, 25, 31, 37))
        _bound(excitation - resisting, terms[:, 43:49], 1e-6,
               f"{prefix} source body{number} reported force sum")


@pytest.mark.parametrize("prefix", ["dense", "seed2_dense",
                                    "step005", "step0025", "step00125",
                                    "step000625", "step0003125"])
def test_dense_source_adjusted_mass_and_joint_force_balance(prefix):
    """Reconstruct source inertia from its body, PTO, and MoorDyn logs."""
    records = [_read(f"{prefix}_body{number}.csv") for number in (1, 2)]
    terms = [_read(f"{prefix}_forces_body{number}.csv") for number in (1, 2)]
    mooring = _read(f"{prefix}_mooring.csv")
    hydro = (Path(APPLICATIONS) /
             "_Common_Input_Files/RM3/hydroData/rm3.h5")
    matrices = []
    with h5py.File(hydro) as h5:
        for number, pitch_inertia in ((1, 21_306_090.66),
                                      (2, 94_407_091.24)):
            group = h5[f"body{number}"]
            full_added = 1000 * np.asarray(
                group["hydro_coeffs/added_mass/inf_freq"])
            added = full_added[:, 6 * (number - 1):6 * number]
            mass = 1000 * float(np.asarray(group["properties/disp_vol"])[0, 0])
            center_z = float(np.asarray(group["properties/cg"]).ravel()[2])
            shift = 2 * np.trace(added[:3, :3])
            adjusted = np.diag([mass + shift] * 3 +
                               [0, pitch_inertia + added[4, 4], 0])
            matrices.append((adjusted, added[4, 4], center_z))

    residuals = np.zeros((1001, 4))
    for index in range(1001):
        angle = records[0][index, 5]
        rate = records[0][index, 11]
        sine, cosine = np.sin(angle), np.cos(angle)
        inertia = np.zeros(4)
        applied = np.zeros(4)
        slider_speed = []
        jacobians = []
        for body_index in (0, 1):
            adjusted, pitch_added, center_z = matrices[body_index]
            pose = records[body_index][index, 1:7]
            velocity = records[body_index][index, 7:13]
            acceleration = terms[body_index][index, 1:7]
            slide = ((pose[2] - center_z
                      - center_z * (cosine - 1)) / cosine)
            radius = center_z + slide
            slider_speed.append((velocity[2] + radius * sine * rate) / cosine)
            jacobian = np.zeros((6, 4))
            jacobian[0, 0] = 1
            jacobian[0, body_index + 1] = sine
            jacobian[0, 3] = radius * cosine
            jacobian[2, body_index + 1] = cosine
            jacobian[2, 3] = -radius * sine
            jacobian[4, 3] = 1
            jacobians.append(jacobian)
            # The source restores rotational, but not translational, mass
            # in its reported forceTotal during postprocessing.
            hydro_force = terms[body_index][index, 43:49].copy()
            hydro_force[4] += pitch_added * acceleration[4]
            inertia += jacobian.T @ adjusted @ acceleration
            applied += jacobian.T @ hydro_force
        pto_force = -1_200_000 * (slider_speed[0] - slider_speed[1])
        applied[1] += pto_force
        applied[2] -= pto_force
        mooring_jacobian = jacobians[1].copy()
        mooring_jacobian[:3, 3] += [21.5 * cosine, 0, -21.5 * sine]
        applied += mooring_jacobian.T @ mooring[index, 13:19]
        residuals[index] = inertia - applied
    _bound(residuals, np.zeros_like(residuals), 1e-3,
           f"{prefix} source adjusted-mass four-coordinate force balance")


def test_saved_source_mooring_load_through_python_body_solver():
    """Separate body integration from MoorDyn state feedback over ten seconds."""
    source_mooring = _read("dense_mooring.csv")

    class PrescribedMoorDyn(MoorDyn):
        def __init__(self, force):
            self.force = force
            self._started = False

        def start(self, position, velocity):
            self._started = True
            return self

        def step(self, position, velocity, time, dt):
            index = round((time + dt) / 0.01)
            return self.force[index].copy()

        def close(self):
            self._started = False

    hydro = (Path(APPLICATIONS) /
             "_Common_Input_Files/RM3/hydroData/rm3.h5")
    wec = WEC("RM3 with saved source mooring load")
    float_body = wec.body("float", hydro, inertia=(0, 21_306_090.66, 0))
    spar = wec.body("spar", hydro, inertia=(0, 94_407_091.24, 0))
    wec.floating_joint(
        float_body, spar, damping=1_200_000,
        moordyn=PrescribedMoorDyn(source_mooring[:, 13:19]),
        moordyn_point=spar.at(0, 0, 21.5),
    )
    result = wec.run(
        JONSWAPWave(2, 8, phase_file=Path(REFERENCE) / "dense_phase.csv",
                    discretization="traditional"),
        dt=0.01, end_time=10, ramp_time=0, radiation_memory=60,
        initial_coordinate={"spar_heave": -0.21},
    )
    assert result.time.shape == (1001,)
    source_wave = _read("dense_wave.csv")
    _bound(result.wave_elevation, source_wave[:, 1], 1e-10,
           "dense diagnostic wave elevation")
    for number, name in ((1, "float"), (2, "spar")):
        saved = _read(f"dense_body{number}.csv")
        assert saved.shape == (1001, 25)
        np.savetxt(Path(REFERENCE) / f"python_saved_force_body{number}.csv",
                   np.column_stack((result.time,
                                    result.bodies[name].position,
                                    result.bodies[name].velocity)),
                   delimiter=",")
        for axis, label in ((0, "surge"), (2, "heave"), (4, "pitch")):
            error = result.bodies[name].position[:, axis] - saved[:, 1 + axis]
            print(f"saved-force {name} {label}: "
                  f"maximum {np.max(np.abs(error)):.8g}, "
                  f"at 0.1 s {error[10]:.8g}, "
                  f"1 s {error[100]:.8g}, at 10 s {error[-1]:.8g}")


def test_live_moordyn_on_same_dense_source_sea(tmp_path):
    """Keep source sea and output grid fixed while restoring mooring feedback."""
    apps = Path(APPLICATIONS)
    source_mooring = _read("dense_mooring.csv")
    hydro = apps / "_Common_Input_Files/RM3/hydroData/rm3.h5"
    input_dir = tmp_path / "Mooring"
    input_dir.mkdir()
    lines = input_dir / "lines.txt"
    shutil.copyfile(apps / "Paraview_Visualization/RM3_MoorDyn_Viz/Mooring/lines.txt",
                    lines)
    wec = WEC("RM3 dense source sea with live MoorDyn")
    float_body = wec.body("float", hydro, inertia=(0, 21_306_090.66, 0))
    spar = wec.body("spar", hydro, inertia=(0, 94_407_091.24, 0))
    wec.floating_joint(
        float_body, spar, damping=1_200_000,
        moordyn=MoorDyn(LIBRARY, lines),
        moordyn_point=spar.at(0, 0, 21.5),
    )
    result = wec.run(
        JONSWAPWave(2, 8, phase_file=Path(REFERENCE) / "dense_phase.csv",
                    discretization="traditional"),
        dt=0.01, end_time=10, ramp_time=0, radiation_memory=60,
        initial_coordinate={"spar_heave": -0.21},
    )
    assert result.time.shape == (1001,)
    _bound(result.wave_elevation, _read("dense_wave.csv")[:, 1],
           1e-10, "live dense-sea wave elevation")
    for number, name in ((1, "float"), (2, "spar")):
        saved = _read(f"dense_body{number}.csv")
        np.savetxt(Path(REFERENCE) / f"python_live_dense_body{number}.csv",
                   np.column_stack((result.time,
                                    result.bodies[name].position,
                                    result.bodies[name].velocity)),
                   delimiter=",")
        for axis, label in ((0, "surge"), (2, "heave"), (4, "pitch")):
            error = result.bodies[name].position[:, axis] - saved[:, 1 + axis]
            print(f"live dense-sea {name} {label}: "
                  f"maximum {np.max(np.abs(error)):.8g}, "
                  f"at 0.1 s {error[10]:.8g}, "
                  f"1 s {error[100]:.8g}, at 10 s {error[-1]:.8g}")
    force = dict(result.raw.extra_outputs)["moordyn_connection_force"]
    np.savetxt(Path(REFERENCE) / "python_live_dense_mooring_force.csv",
               np.column_stack((result.time, force)), delimiter=",")
    for axis, label in ((0, "surge force"), (2, "heave force"),
                        (4, "pitch moment")):
        error = force[:, axis] - source_mooring[:, 13 + axis]
        print(f"live dense-sea mooring {label}: "
              f"maximum {np.max(np.abs(error)):.8g}, "
              f"at 0.1 s {error[10]:.8g}, "
              f"1 s {error[100]:.8g}, at 10 s {error[-1]:.8g}")


def test_physical_fine_step_self_convergence(tmp_path):
    """Check Python step convergence and pair a refined source trajectory."""
    apps = Path(APPLICATIONS)
    hydro = apps / "_Common_Input_Files/RM3/hydroData/rm3.h5"
    line_input = (apps /
                  "Paraview_Visualization/RM3_MoorDyn_Viz/Mooring/lines.txt")
    phase = Path(REFERENCE) / "step0025_phase.csv"

    def run(dt):
        input_dir = tmp_path / f"dt-{dt:g}"
        input_dir.mkdir()
        lines = input_dir / "lines.txt"
        shutil.copyfile(line_input, lines)
        wec = WEC("RM3 MoorDyn physical step convergence")
        float_body = wec.body("float", hydro, inertia=(0, 21_306_090.66, 0))
        spar = wec.body("spar", hydro, inertia=(0, 94_407_091.24, 0))
        wec.floating_joint(
            float_body, spar, damping=1_200_000,
            moordyn=MoorDyn(LIBRARY, lines),
            moordyn_point=spar.at(0, 0, 21.5),
        )
        return wec.run(
            JONSWAPWave(2, 8, phase_file=phase,
                        discretization="traditional"),
            dt=dt, end_time=10, ramp_time=0, radiation_memory=60,
            initial_coordinate={"spar_heave": -0.21},
        )

    coarse = run(0.0025)
    fine = run(0.00125)
    _bound(coarse.time, fine.time[::2], 1e-10,
           "physical solver common times")
    source_wave = _read("step0025_wave.csv")
    for label, result in (("coarse", coarse), ("fine", fine)):
        indices = _indices(source_wave[:, 0], dt=result.time[1], end_time=10)
        _bound(result.wave_elevation[indices], source_wave[:, 1],
               1e-10, f"{label} refined-sea wave")
    for number, name in ((1, "float"), (2, "spar")):
        body_coarse = coarse.bodies[name]
        body_fine = fine.bodies[name]
        _bound(body_coarse.position, body_fine.position[::2], 2.5e-4,
               f"{name} physical position self-convergence")
        _bound(body_coarse.velocity, body_fine.velocity[::2], 1e-4,
               f"{name} physical velocity self-convergence")
        for prefix in ("step0025", "step00125", "step000625",
                       "step0003125"):
            source = _read(f"{prefix}_body{number}.csv")
            for label, result in (("coarse", coarse), ("fine", fine)):
                indices = _indices(source[:, 0], dt=result.time[1],
                                   end_time=10)
                surge_error = (result.bodies[name].position[indices, 0]
                               - source[:, 1])
                print(f"{name} {label} surge against {prefix} source: "
                      f"maximum {np.max(np.abs(surge_error)):.8g} m; "
                      f"at 10 s {surge_error[-1]:.8g} m")
        source = _read(f"step0003125_body{number}.csv")
        indices = _indices(source[:, 0], dt=fine.time[1], end_time=10)
        for axis, label in ((0, "surge"), (2, "heave"), (4, "pitch")):
            _bound(body_fine.position[indices, axis], source[:, 1 + axis],
                   FINE_SOURCE_LIMITS["position"][axis == 4],
                   f"fine-source {name} {label} position")
            _bound(body_fine.velocity[indices, axis], source[:, 7 + axis],
                   FINE_SOURCE_LIMITS["velocity"][axis == 4],
                   f"fine-source {name} {label} velocity")
    source_mooring = _read("step0003125_mooring.csv")
    indices = _indices(source_mooring[:, 0], dt=fine.time[1], end_time=10)
    outputs = dict(fine.raw.extra_outputs)
    for kind, start in (("position", 1), ("velocity", 7), ("force", 13)):
        actual = outputs[f"moordyn_connection_{kind}"][indices]
        for axis, label in ((0, "surge"), (2, "heave"), (4, "pitch")):
            limit = FINE_SOURCE_LIMITS[kind][axis == 4]
            _bound(actual[:, axis], source_mooring[:, start + axis], limit,
                   f"fine-source MoorDyn {label} {kind}")


@pytest.mark.xfail(
    strict=True, raises=PublishedCoupledMotionGap,
    reason="published 0.01 s RM3 MoorDyn coupled motion is not yet paired",
)
def test_public_floating_joint_tracks_published_moordyn_viz(tmp_path):
    apps = Path(APPLICATIONS)
    reference = Path(REFERENCE)
    hydro = apps / "_Common_Input_Files/RM3/hydroData/rm3.h5"
    input_dir = tmp_path / "Mooring"
    input_dir.mkdir()
    lines = input_dir / "lines.txt"
    shutil.copyfile(apps / "Paraview_Visualization/RM3_MoorDyn_Viz/Mooring/lines.txt",
                    lines)

    wec = WEC("Published RM3 MoorDyn visualization")
    float_body = wec.body("float", hydro, inertia=(0, 21_306_090.66, 0))
    spar = wec.body("spar", hydro, inertia=(0, 94_407_091.24, 0))
    wec.floating_joint(
        float_body, spar, damping=1_200_000,
        moordyn=MoorDyn(LIBRARY, lines),
        moordyn_point=spar.at(0, 0, 21.5),
    )
    result = wec.run(
        JONSWAPWave(2, 8, phase_file=reference / "phase.csv",
                    discretization="traditional"),
        dt=0.01, end_time=80, ramp_time=0, radiation_memory=60,
        initial_coordinate={"spar_heave": -0.21},
    )
    assert result.time.shape == (8001,)
    violations = []

    def record(actual, expected, limit, name):
        assert actual.shape == expected.shape, name
        error = float(np.max(np.abs(actual - expected)))
        print(f"{name}: maximum difference {error:.8g}; gate {limit:.8g}")
        if error >= limit:
            violations.append(f"{name}: {error:.8g} >= {limit:.8g}")

    wave = _read("wave.csv")
    _bound(result.wave_elevation[_indices(wave[:, 0])],
           wave[:, 1], 1e-10, "public wave elevation")
    for number, name in ((1, "float"), (2, "spar")):
        saved = _read(f"body{number}.csv")
        indices = _indices(saved[:, 0])
        np.savetxt(reference / f"python_body{number}.csv",
                   np.column_stack((saved[:, 0],
                                    result.bodies[name].position[indices],
                                    result.bodies[name].velocity[indices])),
                   delimiter=",")
        for axis, position_limit, speed_limit in (
            (0, 0.01, 0.005),
            (2, 0.005, 0.005),
            (4, 0.0005, 0.0005),
        ):
            record(result.bodies[name].position[indices, axis],
                   saved[:, 1 + axis], position_limit,
                   f"{name} position axis {axis}")
            record(result.bodies[name].velocity[indices, axis],
                   saved[:, 7 + axis], speed_limit,
                   f"{name} velocity axis {axis}")
    pto = _read("pto.csv")
    pto_force = result.ptos["relative_heave"].force[_indices(pto[:, 0])]
    np.savetxt(reference / "python_pto.csv",
               np.column_stack((pto[:, 0], pto_force)), delimiter=",")
    record(pto_force, pto[:, 15], 2_000, "PTO internal force")
    source_mooring = _read("mooring.csv")
    outputs = dict(result.raw.extra_outputs)
    indices = _indices(source_mooring[:, 0])
    np.savetxt(reference / "python_mooring.csv",
               np.column_stack((source_mooring[:, 0],
                                *(outputs[f"moordyn_connection_{kind}"][indices]
                                  for kind in ("position", "velocity", "force")))),
               delimiter=",")
    for kind, column, limit in (
        ("position", 1, (0.01, 1e-6, 0.005, 1e-6, 0.0005, 1e-6)),
        ("velocity", 7, (0.005, 1e-6, 0.005, 1e-6, 0.0005, 1e-6)),
        ("force", 13, (2_000, 1e-3, 4_000, 1e-3, 2_000, 1e-3)),
    ):
        actual = outputs[f"moordyn_connection_{kind}"][indices]
        for axis, bound in enumerate(limit):
            record(actual[:, axis], source_mooring[:, column + axis], bound,
                   f"mooring {kind} axis {axis}")
    source_tension = _read("fairlead_tension.csv")
    actual_tension = np.loadtxt(input_dir / "lines.out", skiprows=1)
    np.savetxt(reference / "python_fairlead_tension.csv",
               actual_tension[:, :4], delimiter=",")
    record(actual_tension[:, 1:4], source_tension[:, 1:4],
           4_000, "three fairlead tensions")

    # Source Simscape applies added mass through delayed acceleration. Keep
    # this as a diagnostic so the default implicit-mass gate stays physical.
    delay_dir = tmp_path / "MooringDelay"
    delay_dir.mkdir()
    delay_lines = delay_dir / "lines.txt"
    shutil.copyfile(
        apps / "Paraview_Visualization/RM3_MoorDyn_Viz/Mooring/lines.txt",
        delay_lines,
    )
    delayed_wec = WEC("RM3 source added-mass diagnostic")
    delayed_float = delayed_wec.body(
        "float", hydro, inertia=(0, 21_306_090.66, 0),
    )
    delayed_spar = delayed_wec.body(
        "spar", hydro, inertia=(0, 94_407_091.24, 0),
    )
    delayed_wec.floating_joint(
        delayed_float, delayed_spar, damping=1_200_000,
        moordyn=MoorDyn(LIBRARY, delay_lines),
        moordyn_point=delayed_spar.at(0, 0, 21.5),
        added_mass_scheme="simulink_delay",
    )
    try:
        delayed = delayed_wec.run(
            JONSWAPWave(2, 8, phase_file=reference / "phase.csv",
                        discretization="traditional"),
            dt=0.01, end_time=80, ramp_time=0, radiation_memory=60,
            initial_coordinate={"spar_heave": -0.21},
        )
    except (RuntimeError, ValueError) as exc:
        print(f"delayed added-mass diagnostic failed: {exc}")
    else:
        for number, name in ((1, "float"), (2, "spar")):
            saved = _read(f"body{number}.csv")
            indices = _indices(saved[:, 0])
            np.savetxt(reference / f"python_delay_body{number}.csv",
                       np.column_stack((saved[:, 0],
                                        delayed.bodies[name].position[indices],
                                        delayed.bodies[name].velocity[indices])),
                       delimiter=",")
            error = np.max(np.abs(
                delayed.bodies[name].position[indices, 0] - saved[:, 1],
            ))
            print(f"{name} delayed-mass surge position: {error:.8g} m")
        delayed_mooring = dict(delayed.raw.extra_outputs)[
            "moordyn_connection_force"][_indices(source_mooring[:, 0])]
        for axis, label in ((0, "surge force"), (4, "pitch moment")):
            error = np.max(np.abs(
                delayed_mooring[:, axis] - source_mooring[:, 13 + axis],
            ))
            print(f"delayed-mass mooring {label}: {error:.8g}")
    unexpected = {entry.split(": ", 1)[0] for entry in violations}
    unexpected -= EXPECTED_PUBLISHED_GAPS
    assert not unexpected, "; ".join(violations)
    if violations:
        raise PublishedCoupledMotionGap("; ".join(violations))

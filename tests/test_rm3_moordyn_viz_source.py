"""Pair the published Traditional JONSWAP RM3 MoorDyn visualization motion."""

import os
from pathlib import Path
import shutil

import numpy as np
import pytest

from wecsim import JONSWAPWave, MoorDyn, WEC
from wecsim.irregularWave import (
    jonswap_equal_energy_components, synthesize_irregular_response,
)


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


def _indices(source_time, *, dt=0.01, end_time=80):
    indices = np.rint(source_time / dt).astype(int)
    assert np.all((indices >= 0) & (indices <= round(end_time / dt)))
    np.testing.assert_allclose(indices * dt, source_time, rtol=0, atol=1e-8)
    return indices


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
        for axis, label in ((0, "surge"), (2, "heave"), (4, "pitch")):
            error = dense[:, 1 + axis] - coarse[:, 1 + axis]
            print(f"source output-step {name} {label}: "
                  f"maximum {np.max(np.abs(error)):.8g}, "
                  f"at 1 s {error[10]:.8g}, at 10 s {error[-1]:.8g}")
    dense_mooring = _read("dense_mooring.csv")[::10]
    coarse_mooring = _read("coarse_mooring.csv")
    assert dense_mooring.shape == coarse_mooring.shape == (101, 19)
    for axis, label in ((0, "surge force"), (2, "heave force"),
                        (4, "pitch moment")):
        error = dense_mooring[:, 13 + axis] - coarse_mooring[:, 13 + axis]
        print(f"source output-step mooring {label}: "
              f"maximum {np.max(np.abs(error)):.8g}, "
              f"at 1 s {error[10]:.8g}, at 10 s {error[-1]:.8g}")


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
    record(result.wave_elevation[_indices(wave[:, 0])],
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
    assert not violations, "; ".join(violations)

"""Isolate source solver maximum step from MoorDyn coupling step."""

import os
from pathlib import Path
import shutil
from zipfile import ZipFile

import numpy as np
import pytest

from wecsim import JONSWAPWave, MoorDyn, WEC


RECORDS = os.environ.get("WEC_SIM_RM3_SOLVER_COUPLING_AUDIT_DIR")
SOURCE = os.environ.get("WEC_SIM_MATLAB_SOURCE_DIR")
APPLICATIONS = os.environ.get("WEC_SIM_APPLICATIONS_DIR")
LIBRARY = os.environ.get("WEC_SIM_MOORDYN_LIBRARY")
pytestmark = pytest.mark.skipif(
    not (RECORDS and SOURCE), reason="pinned solver/coupling audit absent",
)


def _read(label, channel):
    values = np.loadtxt(Path(RECORDS) / f"{label}_{channel}.csv",
                        delimiter=",", ndmin=2)
    assert np.isfinite(values).all(), f"{label} {channel}"
    return values


def test_source_moordyn_coupling_pulse_uses_simu_dt():
    model = (Path(SOURCE) / "lib/WEC-Sim/WECSim_Lib_Moorings.slx")
    with ZipFile(model) as archive:
        subsystem = archive.read(
            "simulink/systems/system_12_49.xml",
        ).decode()
    assert '<P Name="Value">simu.dt</P>' in subsystem
    assert '<P Name="Period">2*simu.dt</P>' in subsystem
    assert '<P Name="PulseWidth">50</P>' in subsystem


def test_independent_solver_and_moordyn_step_audit():
    labels = ("published_step", "solver_only", "solver_and_moordyn")
    channels = ("body1", "body2", "mooring")
    phase = [_read(label, "phase") for label in labels]
    wave = [_read(label, "wave") for label in labels]
    for values in phase:
        assert values.shape == (1000, 1)
        np.testing.assert_allclose(values, phase[0], rtol=0, atol=1e-13)
    for values in wave:
        assert values.shape == (1001, 2)
        np.testing.assert_allclose(values[:, 0], np.arange(1001) * .01,
                                   rtol=0, atol=1e-9)
        np.testing.assert_allclose(values, wave[0], rtol=0, atol=1e-10)

    for channel in channels:
        records = [_read(label, channel) for label in labels]
        shape = (1001, 19) if channel == "mooring" else (1001, 13)
        assert all(record.shape == shape for record in records)
        for record in records:
            np.testing.assert_allclose(record[:, 0], wave[0][:, 0],
                                       rtol=0, atol=1e-9)
            np.testing.assert_allclose(record[0, 1:13], records[0][0, 1:13],
                                       rtol=0, atol=1e-9)
        # Quantify which independent refinement removes the source drift.
        # No outcome is pre-assumed; the observed numerical bounds are
        # recorded after the pinned source experiment runs.
        for name, index in (("surge", 1), ("heave", 3), ("pitch", 5)):
            solver_change = np.max(abs(records[1][:, index] - records[0][:, index]))
            dt_change = np.max(abs(records[2][:, index] - records[1][:, index]))
            print(f"{channel} {name}: MaxStep-only {solver_change:.9g}; "
                  f"additional simu.dt change {dt_change:.9g}")
            if channel in ("body1", "body2") and name == "surge":
                # Pinned source result: the MaxStep-only arm removes nearly
                # all of the published-step surge shift. The remaining
                # simu.dt change also refines wave/radiation grids, so it
                # must not be attributed solely to MoorDyn coupling.
                assert solver_change > .17
                assert dt_change < .001
            if channel in ("body1", "body2") and name == "pitch":
                assert solver_change > .0017
                assert dt_change < 2e-5
        if channel == "mooring":
            for name, index in (("surge force", 13), ("heave force", 15),
                                ("pitch moment", 17)):
                solver_change = np.max(abs(records[1][:, index] - records[0][:, index]))
                dt_change = np.max(abs(records[2][:, index] - records[1][:, index]))
                print(f"MoorDyn {name}: MaxStep-only {solver_change:.9g}; "
                      f"additional simu.dt change {dt_change:.9g}")


@pytest.mark.skipif(not (RECORDS and APPLICATIONS and LIBRARY),
                    reason="pinned source, hydro data, or MoorDyn absent")
def test_independent_python_motion_on_solver_refined_source(tmp_path):
    """Advance Python with source sea phases, never source forces or motion."""
    apps = Path(APPLICATIONS)
    hydro = apps / "_Common_Input_Files/RM3/hydroData/rm3.h5"
    lines = tmp_path / "Mooring" / "lines.txt"
    lines.parent.mkdir()
    shutil.copyfile(apps / "Paraview_Visualization/RM3_MoorDyn_Viz"
                    / "Mooring/lines.txt", lines)
    wec = WEC("RM3 solver-only MATLAB step pair")
    float_body = wec.body("float", hydro, inertia=(0, 21_306_090.66, 0))
    spar = wec.body("spar", hydro, inertia=(0, 94_407_091.24, 0))
    wec.floating_joint(
        float_body, spar, damping=1_200_000,
        moordyn=MoorDyn(LIBRARY, lines),
        moordyn_point=spar.at(0, 0, 21.5),
    )
    result = wec.run(
        JONSWAPWave(
            2, 8, phase_file=Path(RECORDS) / "solver_only_phase.csv",
            discretization="traditional",
        ),
        dt=.00125, end_time=10, ramp_time=0, radiation_memory=60,
        initial_coordinate={"spar_heave": -.21},
    )
    indices = np.arange(1001) * 8
    np.testing.assert_allclose(result.time[indices], _read("solver_only", "wave")[:, 0],
                               rtol=0, atol=1e-8)
    np.testing.assert_allclose(result.wave_elevation[indices],
                               _read("solver_only", "wave")[:, 1],
                               rtol=0, atol=1e-10)

    limits = {"position": (0.002, 5e-5),
              "velocity": (0.0007, 5e-5),
              "force": (1500, 10000)}
    for number, name in ((1, "float"), (2, "spar")):
        source = _read("solver_only", f"body{number}")
        for kind, offset in (("position", 1), ("velocity", 7)):
            actual = getattr(result.bodies[name], kind)[indices]
            for axis in (0, 2, 4):
                error = np.max(abs(actual[:, axis] - source[:, offset + axis]))
                assert error < limits[kind][axis == 4], (name, kind, axis, error)

    source = _read("solver_only", "mooring")
    outputs = dict(result.raw.extra_outputs)
    for kind, offset in (("position", 1), ("velocity", 7), ("force", 13)):
        actual = outputs[f"moordyn_connection_{kind}"][indices]
        for axis in (0, 2, 4):
            error = np.max(abs(actual[:, axis] - source[:, offset + axis]))
            assert error < limits[kind][axis == 4], ("MoorDyn", kind, axis, error)

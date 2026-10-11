"""Isolate source solver maximum step from MoorDyn coupling step."""

import os
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pytest


RECORDS = os.environ.get("WEC_SIM_RM3_SOLVER_COUPLING_AUDIT_DIR")
SOURCE = os.environ.get("WEC_SIM_MATLAB_SOURCE_DIR")
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
            coupling_change = np.max(abs(records[2][:, index] - records[1][:, index]))
            print(f"{channel} {name}: MaxStep-only {solver_change:.9g}; "
                  f"additional MoorDyn coupling {coupling_change:.9g}")
        if channel == "mooring":
            for name, index in (("surge force", 13), ("heave force", 15),
                                ("pitch moment", 17)):
                solver_change = np.max(abs(records[1][:, index] - records[0][:, index]))
                coupling_change = np.max(abs(records[2][:, index] - records[1][:, index]))
                print(f"MoorDyn {name}: MaxStep-only {solver_change:.9g}; "
                      f"additional coupling {coupling_change:.9g}")

"""Diagnose full-duration RM3 motion with MATLAB ode45 MaxStep refined alone."""

import os
from pathlib import Path
import shutil

import numpy as np
import pytest

from wecsim import JONSWAPWave, MoorDyn, WEC


APPLICATIONS = os.environ.get("WEC_SIM_APPLICATIONS_DIR")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_RM3_FULL_SOLVER_STEP_DIR")
LIBRARY = os.environ.get("WEC_SIM_MOORDYN_LIBRARY")
pytestmark = pytest.mark.skipif(
    not (APPLICATIONS and REFERENCE and LIBRARY),
    reason="full-duration MaxStep-only MATLAB run and pinned MoorDyn absent",
)


def _load(name):
    values = np.loadtxt(Path(REFERENCE) / f"{name}.csv", delimiter=",", ndmin=2)
    assert np.isfinite(values).all(), name
    return values


def _difference(actual, expected, label):
    assert actual.shape == expected.shape, label
    error = float(np.max(np.abs(actual - expected)))
    print(f"{label}: {error:.9g}")
    return error


def test_full_duration_maxstep_only_source_and_python(tmp_path):
    phase = _load("phase")
    wave = _load("wave")
    assert phase.shape == (1000, 1)
    assert wave.shape == (8001, 2)
    _difference(wave[:, 0], np.arange(8001) * .01, "source output time")
    assert np.max(abs(wave[:, 0] - np.arange(8001) * .01)) < 1e-8

    apps = Path(APPLICATIONS)
    hydro = apps / "_Common_Input_Files/RM3/hydroData/rm3.h5"
    lines = tmp_path / "Mooring" / "lines.txt"
    lines.parent.mkdir()
    shutil.copyfile(apps / "Paraview_Visualization/RM3_MoorDyn_Viz"
                    / "Mooring/lines.txt", lines)
    wec = WEC("RM3 80 s source MaxStep-only comparison")
    float_body = wec.body("float", hydro, inertia=(0, 21_306_090.66, 0))
    spar = wec.body("spar", hydro, inertia=(0, 94_407_091.24, 0))
    wec.floating_joint(
        float_body, spar, damping=1_200_000,
        moordyn=MoorDyn(LIBRARY, lines),
        moordyn_point=spar.at(0, 0, 21.5),
    )
    result = wec.run(
        JONSWAPWave(2, 8, phase_file=Path(REFERENCE) / "phase.csv",
                    discretization="traditional"),
        dt=.00125, end_time=80, ramp_time=0, radiation_memory=60,
        initial_coordinate={"spar_heave": -.21},
    )
    indices = np.arange(8001) * 8
    assert result.time.shape == (64001,)
    assert _difference(result.time[indices], wave[:, 0],
                       "Python sample time") < 1e-8
    assert _difference(result.wave_elevation[indices], wave[:, 1],
                       "Python wave elevation") < 1e-10
    for number, name in ((1, "float"), (2, "spar")):
        source = _load(f"body{number}")
        assert source.shape == (8001, 25)
        for kind, offset in (("position", 1), ("velocity", 7)):
            actual = getattr(result.bodies[name], kind)[indices]
            for axis, label in ((0, "surge"), (2, "heave"), (4, "pitch")):
                _difference(actual[:, axis], source[:, offset + axis],
                            f"{name} {label} {kind}")
    source = _load("mooring")
    assert source.shape == (8001, 19)
    outputs = dict(result.raw.extra_outputs)
    for kind, offset in (("position", 1), ("velocity", 7), ("force", 13)):
        actual = outputs[f"moordyn_connection_{kind}"][indices]
        for axis, label in ((0, "surge"), (2, "heave"), (4, "pitch")):
            _difference(actual[:, axis], source[:, offset + axis],
                        f"MoorDyn {label} {kind}")
    source = _load("pto")
    assert source.shape == (8001, 25)
    _difference(result.ptos["relative_heave"].force[indices], source[:, 15],
                "PTO force")

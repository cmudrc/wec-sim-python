"""Full-duration fine-step RM3 MoorDyn source and independent motion gate."""

import os
from pathlib import Path
import shutil

import numpy as np
import pytest

from wecsim import JONSWAPWave, MoorDyn, WEC


APPLICATIONS = os.environ.get("WEC_SIM_APPLICATIONS_DIR")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_RM3_VIZ_FULL_FINE_DIR")
LIBRARY = os.environ.get("WEC_SIM_MOORDYN_LIBRARY")
pytestmark = pytest.mark.skipif(
    not (APPLICATIONS and REFERENCE and LIBRARY),
    reason="full-duration fine-step MATLAB output and pinned MoorDyn absent",
)

LIMITS = {
    "position": (0.002, 5e-5),
    "velocity": (0.0005, 5e-5),
    "force": (500, 2000),
}


def _read(label, kind):
    values = np.loadtxt(Path(REFERENCE) / f"{label}_{kind}.csv",
                        delimiter=",", ndmin=2)
    assert np.isfinite(values).all(), f"{label} {kind}"
    return values


def _compare(actual, expected, limit, label, failures):
    assert actual.shape == expected.shape, label
    error = float(np.max(np.abs(actual - expected)))
    print(f"{label}: {error:.8g} (gate {limit:.8g})")
    if not np.isfinite(error) or error >= limit:
        failures.append(f"{label}: {error:.8g} >= {limit:.8g}")


def _active_channels(actual, expected, kinds, label, failures):
    for kind, first, axes in kinds:
        for axis, name in axes:
            _compare(actual[:, first + axis], expected[:, first + axis],
                     LIMITS[kind][axis == 4],
                     f"{label} {name} {kind}", failures)


BODY_KINDS = (
    ("position", 1, ((0, "surge"), (2, "heave"), (4, "pitch"))),
    ("velocity", 7, ((0, "surge"), (2, "heave"), (4, "pitch"))),
)
MOORING_KINDS = BODY_KINDS + (
    ("force", 13, ((0, "surge"), (2, "heave"), (4, "pitch"))),
)


def test_full_duration_refined_source_and_python_motion(tmp_path):
    coarse, fine = "fullfine_000625", "fullfine_0003125"
    failures = []
    coarse_phase = _read(coarse, "phase")
    fine_phase = _read(fine, "phase")
    assert coarse_phase.shape == fine_phase.shape == (1000, 1)
    _compare(coarse_phase, fine_phase, 1e-13, "source phase", failures)
    coarse_wave, fine_wave = _read(coarse, "wave"), _read(fine, "wave")
    assert coarse_wave.shape == fine_wave.shape == (8001, 2)
    _compare(coarse_wave, fine_wave, 1e-10, "source wave", failures)
    _compare(fine_wave[:, 0], np.arange(8001) * .01, 1e-8,
             "source sample time", failures)

    for number, name in ((1, "float"), (2, "spar")):
        left, right = (_read(prefix, f"body{number}")
                       for prefix in (coarse, fine))
        assert left.shape == right.shape == (8001, 25)
        _active_channels(left, right, BODY_KINDS,
                         f"source step change {name}", failures)
    left, right = (_read(prefix, "mooring")
                   for prefix in (coarse, fine))
    assert left.shape == right.shape == (8001, 19)
    _active_channels(left, right, MOORING_KINDS,
                     "source step change MoorDyn", failures)
    left, right = (_read(prefix, "pto") for prefix in (coarse, fine))
    assert left.shape == right.shape == (8001, 25)
    _compare(left[:, 15], right[:, 15], 2000,
             "source step change PTO force", failures)

    apps = Path(APPLICATIONS)
    hydro = apps / "_Common_Input_Files/RM3/hydroData/rm3.h5"
    lines = tmp_path / "Mooring" / "lines.txt"
    lines.parent.mkdir()
    shutil.copyfile(apps / "Paraview_Visualization/RM3_MoorDyn_Viz"
                    / "Mooring/lines.txt", lines)
    wec = WEC("RM3 MoorDyn 80 s fine-step comparison")
    float_body = wec.body("float", hydro, inertia=(0, 21_306_090.66, 0))
    spar = wec.body("spar", hydro, inertia=(0, 94_407_091.24, 0))
    wec.floating_joint(
        float_body, spar, damping=1_200_000,
        moordyn=MoorDyn(LIBRARY, lines),
        moordyn_point=spar.at(0, 0, 21.5),
    )
    result = wec.run(
        JONSWAPWave(2, 8, phase_file=Path(REFERENCE) / f"{fine}_phase.csv",
                    discretization="traditional"),
        dt=.00125, end_time=80, ramp_time=0, radiation_memory=60,
        initial_coordinate={"spar_heave": -.21},
    )
    assert result.time.shape == (64001,)
    indices = np.arange(8001) * 8
    _compare(result.time[indices], fine_wave[:, 0], 1e-8,
             "Python sample time", failures)
    _compare(result.wave_elevation[indices], fine_wave[:, 1], 1e-10,
             "Python wave elevation", failures)
    for number, name in ((1, "float"), (2, "spar")):
        source = _read(fine, f"body{number}")
        actual = np.column_stack((
            result.time[indices], result.bodies[name].position[indices],
            result.bodies[name].velocity[indices],
        ))
        _active_channels(actual, source, BODY_KINDS,
                         f"Python versus fine source {name}", failures)
    source = _read(fine, "mooring")
    outputs = dict(result.raw.extra_outputs)
    actual = np.column_stack((
        result.time[indices],
        *(outputs[f"moordyn_connection_{kind}"][indices]
          for kind in ("position", "velocity", "force")),
    ))
    _active_channels(actual, source, MOORING_KINDS,
                     "Python versus fine source MoorDyn", failures)
    source = _read(fine, "pto")
    _compare(result.ptos["relative_heave"].force[indices], source[:, 15],
             2000, "Python versus fine source PTO force", failures)
    assert not failures, "; ".join(failures)

"""Gate the stable 120 s prefix and audit a longer selected-heading path."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim import PMWave, WEC, WorldPoint


APPLICATIONS = os.environ.get("WEC_SIM_APPLICATIONS_DIR")
SOURCE = os.environ.get("WEC_SIM_MATLAB_250S_AUDIT_DIR")
pytestmark = pytest.mark.skipif(
    not (APPLICATIONS and SOURCE),
    reason="pinned variable-yaw quarter-step source output not provided",
)


def _errors(fine, quarter, result):
    """Compare a quarter-step prefix with the same half-step source window."""
    fine = fine[:(len(quarter) + 1) // 2]
    motion = result.bodies["flap"]
    pto = result.ptos["hinge"]
    source_work = -np.trapezoid(quarter[:, 6], quarter[:, 0])
    fine_work = -np.trapezoid(fine[:, 6], fine[:, 0])
    python_work = np.trapezoid(pto.absorbed_power[:len(quarter)],
                               result.time[:len(quarter)])
    assert source_work > 0 and fine_work > 0 and python_work > 0
    source_step = (
        np.max(np.abs(fine[:, 1] - quarter[::2, 1])),
        np.max(np.abs(fine[:, 2] - quarter[::2, 2])),
        np.max(np.abs(fine[:, 5] - quarter[::2, 5])),
        abs(source_work - fine_work) / source_work,
    )
    python_source = (
        np.max(np.abs(motion.position[:len(quarter), 5] - quarter[:, 1])),
        np.max(np.abs(motion.velocity[:len(quarter), 5] - quarter[:, 2])),
        np.max(np.abs(pto.force[:len(quarter)] - quarter[:, 5])),
        abs(python_work - source_work) / source_work,
    )
    return source_step, python_source


def _report(duration, source_step, python_source):
    print(
        f"{duration} s quarter-step source self-change: "
        f"yaw={source_step[0]:.6g} rad, speed={source_step[1]:.6g} rad/s, "
        f"torque={source_step[2]:.6g} Nm, work={source_step[3]:.3%}; "
        f"Python/source yaw={python_source[0]:.6g} rad, "
        f"speed={python_source[1]:.6g} rad/s, "
        f"torque={python_source[2]:.6g} Nm, work={python_source[3]:.3%}"
    )


def test_irregular_140s_quarter_step_audit():
    source = Path(SOURCE)
    quarter_file = source / "quarter_trajectory.csv"
    assert quarter_file.is_file(), "quarter-step source run is absent"
    fine = np.loadtxt(source / "fine_trajectory.csv", delimiter=",")
    quarter = np.loadtxt(quarter_file, delimiter=",")
    fine_wave = np.loadtxt(source / "fine_wave.csv", delimiter=",")
    quarter_wave = np.loadtxt(source / "quarter_wave.csv", delimiter=",")
    assert fine.shape == (50001, 7)
    assert quarter.shape == (56001, 7)
    assert fine_wave.shape == (50001, 2)
    assert quarter_wave.shape == (56001, 2)
    assert all(np.isfinite(array).all() for array in
               (fine, quarter, fine_wave, quarter_wave))
    fine140 = fine[:28001]
    fine140_wave = fine_wave[:28001]
    np.testing.assert_allclose(quarter[::2, 0], fine140[:, 0],
                               rtol=0, atol=1e-10)
    np.testing.assert_allclose(quarter_wave[::2], fine140_wave, rtol=0,
                               atol=1e-12)

    headings = np.arange(-30, 30.25, .25)
    relative = 10 - np.degrees(quarter[:, 1])
    selected = headings[np.abs(relative[:, None] - headings).argmin(axis=1)]
    np.testing.assert_array_equal(quarter[:, 4], selected)
    np.testing.assert_allclose(quarter[:, 5], -120000 * quarter[:, 2],
                               rtol=0, atol=1e-8)
    np.testing.assert_allclose(quarter[:, 6], quarter[:, 5] * quarter[:, 2],
                               rtol=0, atol=1e-8)

    hydro = (Path(APPLICATIONS) /
             "_Common_Input_Files/OSWEC/hydroData/oswec.h5")
    wec = WEC("OSWEC 140 s selected-heading diagnostic")
    flap = wec.body("flap", hydro, mass=12700, inertia=(1.85e6,) * 3,
                    passive_yaw=True, yaw_heading_bank=headings)
    wec.body("base", hydro, mass=999, inertia=(999,) * 3)
    yaw = wec.coordinate("yaw", flap.move("yaw", pivot=WorldPoint(0, 0, -8.9)))
    wec.rotational_pto("hinge", yaw, damping=120000)
    result = wec.run(
        PMWave(2.5, 8, direction=10, seed=1, phase_generator="matlab"),
        dt=.0025, end_time=140, ramp_time=100, radiation_memory=40,
    )
    np.testing.assert_allclose(result.time, quarter[:, 0], rtol=0, atol=1e-10)
    np.testing.assert_allclose(result.wave_elevation, quarter_wave[:, 1],
                               rtol=0, atol=1e-12)

    # These limits exceed the measured 120 s 0.01-to-0.005 s source
    # self-change but still reject a new event-driven trajectory split.
    limits = (.0006, .00015, 20, .001)
    source120, python120 = _errors(fine, quarter[:48001], result)
    _report(120, source120, python120)
    for measured in (source120, python120):
        assert all(value < limit for value, limit in zip(measured, limits))

    # The pinned 0.01-to-0.005 s source pair stays within these same limits
    # through 140 s, then exceeds the yaw limit by 145 s. The quarter-step
    # pair must establish its own stability before extending the gate.
    source140, python140 = _errors(fine, quarter, result)
    _report(140, source140, python140)
    source_stable = all(value < limit for value, limit in zip(source140, limits))
    print(f"140 s MATLAB step stability under 120 s limits: {source_stable}")
    assert source_stable, "140 s MATLAB source is step-sensitive"
    assert all(value < limit for value, limit in zip(python140, limits))

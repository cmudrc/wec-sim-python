"""Gate the early selected-heading motion against a refined source path."""

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


def test_irregular_120s_quarter_step_pair():
    source = Path(SOURCE)
    quarter_file = source / "quarter_trajectory.csv"
    assert quarter_file.is_file(), "quarter-step source run is absent"
    fine = np.loadtxt(source / "fine_trajectory.csv", delimiter=",")
    quarter = np.loadtxt(quarter_file, delimiter=",")
    fine_wave = np.loadtxt(source / "fine_wave.csv", delimiter=",")
    quarter_wave = np.loadtxt(source / "quarter_wave.csv", delimiter=",")
    assert fine.shape == (50001, 7)
    assert quarter.shape == (48001, 7)
    assert fine_wave.shape == (50001, 2)
    assert quarter_wave.shape == (48001, 2)
    assert all(np.isfinite(array).all() for array in
               (fine, quarter, fine_wave, quarter_wave))
    fine120 = fine[:24001]
    fine120_wave = fine_wave[:24001]
    np.testing.assert_allclose(quarter[::2, 0], fine120[:, 0], rtol=0, atol=1e-10)
    np.testing.assert_allclose(quarter_wave[::2], fine120_wave, rtol=0,
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
    wec = WEC("OSWEC quarter-step selected-heading diagnostic")
    flap = wec.body("flap", hydro, mass=12700, inertia=(1.85e6,) * 3,
                    passive_yaw=True, yaw_heading_bank=headings)
    wec.body("base", hydro, mass=999, inertia=(999,) * 3)
    yaw = wec.coordinate("yaw", flap.move("yaw", pivot=WorldPoint(0, 0, -8.9)))
    wec.rotational_pto("hinge", yaw, damping=120000)
    result = wec.run(
        PMWave(2.5, 8, direction=10, seed=1, phase_generator="matlab"),
        dt=.0025, end_time=120, ramp_time=100, radiation_memory=40,
    )
    np.testing.assert_allclose(result.time, quarter[:, 0], rtol=0, atol=1e-10)
    np.testing.assert_allclose(result.wave_elevation, quarter_wave[:, 1],
                               rtol=0, atol=1e-12)

    source_work = -np.trapezoid(quarter[:, 6], quarter[:, 0])
    python_work = np.trapezoid(result.ptos["hinge"].absorbed_power, result.time)
    assert source_work > 0 and python_work > 0
    source_step_yaw = np.max(np.abs(fine120[:, 1] - quarter[::2, 1]))
    source_step_speed = np.max(np.abs(fine120[:, 2] - quarter[::2, 2]))
    source_step_torque = np.max(np.abs(fine120[:, 5] - quarter[::2, 5]))
    source_step_work = -np.trapezoid(fine120[:, 6], fine120[:, 0])
    motion = result.bodies["flap"]
    python_source_yaw = np.max(np.abs(motion.position[:, 5] - quarter[:, 1]))
    python_source_speed = np.max(np.abs(motion.velocity[:, 5] - quarter[:, 2]))
    python_source_torque = np.max(np.abs(
        result.ptos["hinge"].force - quarter[:, 5]
    ))
    source_step_work_error = abs(source_work - source_step_work) / source_work
    python_work_error = abs(python_work - source_work) / source_work

    # These limits exceed the measured 120 s 0.01-to-0.005 s source
    # self-change but still reject a new event-driven trajectory split.
    limits = (.0006, .00015, 20, .001)
    print(
        f"120 s quarter-step source self-change: yaw={source_step_yaw:.6g} rad, "
        f"speed={source_step_speed:.6g} rad/s, "
        f"torque={source_step_torque:.6g} Nm, work={source_step_work_error:.3%}; "
        f"Python/source yaw={python_source_yaw:.6g} rad, "
        f"speed={python_source_speed:.6g} rad/s, "
        f"torque={python_source_torque:.6g} Nm, work={python_work_error:.3%}"
    )
    for measured in (
        (source_step_yaw, source_step_speed, source_step_torque,
         source_step_work_error),
        (python_source_yaw, python_source_speed, python_source_torque,
         python_work_error),
    ):
        assert all(value < limit for value, limit in zip(measured, limits))

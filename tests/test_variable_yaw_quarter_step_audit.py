"""Audit selected-heading step sensitivity before claiming long-run parity."""

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


def test_irregular_250s_quarter_step_diagnostic():
    source = Path(SOURCE)
    quarter_file = source / "quarter_trajectory.csv"
    assert quarter_file.is_file(), "quarter-step source run is absent"
    fine = np.loadtxt(source / "fine_trajectory.csv", delimiter=",")
    quarter = np.loadtxt(quarter_file, delimiter=",")
    fine_wave = np.loadtxt(source / "fine_wave.csv", delimiter=",")
    quarter_wave = np.loadtxt(source / "quarter_wave.csv", delimiter=",")
    assert fine.shape == (50001, 7)
    assert quarter.shape == (100001, 7)
    assert fine_wave.shape == (50001, 2)
    assert quarter_wave.shape == (100001, 2)
    assert all(np.isfinite(array).all() for array in
               (fine, quarter, fine_wave, quarter_wave))
    np.testing.assert_allclose(quarter[::2, 0], fine[:, 0], rtol=0, atol=1e-10)
    np.testing.assert_allclose(quarter_wave[::2], fine_wave, rtol=0,
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
        dt=.0025, end_time=250, ramp_time=100, radiation_memory=40,
    )
    np.testing.assert_allclose(result.time, quarter[:, 0], rtol=0, atol=1e-10)
    np.testing.assert_allclose(result.wave_elevation, quarter_wave[:, 1],
                               rtol=0, atol=1e-12)

    source_work = -np.trapezoid(quarter[:, 6], quarter[:, 0])
    python_work = np.trapezoid(result.ptos["hinge"].absorbed_power, result.time)
    assert source_work > 0 and python_work > 0
    source_step_yaw = np.max(np.abs(fine[:, 1] - quarter[::2, 1]))
    source_step_work = -np.trapezoid(fine[:, 6], fine[:, 0])
    python_source_yaw = np.max(np.abs(
        result.bodies["flap"].position[:, 5] - quarter[:, 1]
    ))
    # The discontinuous heading bank can amplify one-sample switch timing.
    # Record source self-change and cross-language difference separately;
    # neither is a trajectory parity gate until source convergence is shown.
    print(
        f"250 s quarter-step source self-change: yaw={source_step_yaw:.6g} rad, "
        f"work={abs(source_work - source_step_work) / source_work:.3%}; "
        f"Python/source yaw={python_source_yaw:.6g} rad, "
        f"work={abs(python_work - source_work) / source_work:.3%}"
    )

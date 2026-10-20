"""Pair hinged flap body damping with a PM sea and radiation memory."""

from pathlib import Path

import numpy as np

from wecsim import PMWave, WEC, WorldPoint
from wecsim.environment.irregularWave import pm_equal_energy_components


ROOT = Path(__file__).resolve().parents[1]
HYDRO = ROOT / "tests/fixtures/matlab/oswec_multi_wave/oswec.h5"
SOURCE = ROOT / "tests/fixtures/matlab/oswec_hinge_body_damping_pm/source.npz"
AXES = (0, 2, 4)
# Retain the predeclared regular-wave limits for the same damping matrix.
POSITION_LIMITS = (0.03, 0.012, 0.007)
SPEED_LIMITS = (0.025, 0.01, 0.006)
DELTA_POSITION_LIMITS = (0.001, 0.0005, 0.00015)
DELTA_SPEED_LIMITS = (0.001, 0.0005, 0.00025)


def _maximum(left, right):
    assert left.shape == right.shape
    return float(np.max(np.abs(left - right)))


def _run(matrix):
    wec = WEC("OSWEC PM hinged body damping")
    flap = wec.body("flap", HYDRO, mass=127_000, inertia=(0, 1.85e6, 0),
                    linear_damping=matrix)
    base = wec.fixed_body("base", center_gravity=(0, 0, -10.9),
                          mass=999, inertia=(1, 1, 1))
    wec.fixed_hinge(flap, base, location=WorldPoint(0, 0, -10),
                    pto_location=WorldPoint(0, 0, -8.9))
    return wec.run(
        PMWave(2.5, 8, direction=0, seed=1, phase_generator="matlab"),
        dt=.025, end_time=150, ramp_time=100, radiation_memory=40,
    )


def test_pm_body_damping_matches_source_change_with_memory_radiation():
    with np.load(SOURCE, allow_pickle=False) as saved:
        flap = saved["body1"].copy()
        base = saved["body2"].copy()
        pto = saved["pto"].copy()
        wave = saved["wave"].copy()
        phase = saved["phase"].copy()
        damping = saved["damping"].copy()
    assert flap.shape == base.shape == (2, 6001, 13)
    assert pto.shape == (2, 6001, 19)
    assert wave.shape == (2, 6001, 2)
    assert phase.shape == (2, 500)
    np.testing.assert_array_equal(damping[0], np.zeros((6, 6)))
    np.testing.assert_array_equal(
        damping[1], np.diag([100_000, 0, 0, 0, 2_000_000, 0]),
    )
    np.testing.assert_array_equal(phase[0], phase[1])
    np.testing.assert_array_equal(wave[0], wave[1])
    assert _maximum(flap[1, :, 5], flap[0, :, 5]) > .02
    assert np.max(np.abs(pto[:, :, 17])) < 1e-6

    components = pm_equal_energy_components(
        HYDRO, significant_height=2.5, peak_period=8,
        directions=[0], spreading=[1], seed=1, phase_generator="matlab",
    )
    np.testing.assert_allclose(components.phase[:, 0], phase[0],
                               rtol=0, atol=1e-12)
    control, damped = _run(None), _run(damping[1])
    for index, result in enumerate((control, damped)):
        np.testing.assert_allclose(result.time, flap[index, :, 0],
                                   rtol=0, atol=1e-10)
        assert _maximum(result.wave_elevation, wave[index, :, 1]) < 1e-12
        actual = result.bodies["flap"]
        for axis, position_limit, speed_limit in zip(
                AXES, POSITION_LIMITS, SPEED_LIMITS, strict=True):
            assert _maximum(actual.position[:, axis],
                            flap[index, :, 1 + axis]) < position_limit
            assert _maximum(actual.velocity[:, axis],
                            flap[index, :, 7 + axis]) < speed_limit
        assert _maximum(result.bodies["base"].position,
                        base[index, :, 1:7]) < 1e-10
        assert _maximum(result.bodies["base"].velocity,
                        base[index, :, 7:13]) < 1e-10
        assert _maximum(result.ptos["hinge"].force,
                        pto[index, :, 17]) < 1e-6

    # Compare the effect of damping to cancel shared solver error. The source
    # arms have identical waves, geometry, and PTO settings.
    for axis, position_limit, speed_limit in zip(
            AXES, DELTA_POSITION_LIMITS, DELTA_SPEED_LIMITS, strict=True):
        source_position_change = flap[1, :, 1 + axis] - flap[0, :, 1 + axis]
        python_position_change = (damped.bodies["flap"].position[:, axis]
                                  - control.bodies["flap"].position[:, axis])
        assert _maximum(python_position_change,
                        source_position_change) < position_limit
        source_speed_change = flap[1, :, 7 + axis] - flap[0, :, 7 + axis]
        python_speed_change = (damped.bodies["flap"].velocity[:, axis]
                               - control.bodies["flap"].velocity[:, axis])
        assert _maximum(python_speed_change,
                        source_speed_change) < speed_limit

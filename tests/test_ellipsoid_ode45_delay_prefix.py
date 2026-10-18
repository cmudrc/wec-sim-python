"""Pin the active added-mass feedback in the published ode45 ellipsoid."""

from pathlib import Path

import numpy as np


FIXTURE = (Path(__file__).parent / "fixtures/matlab"
           / "ellipsoid_ode45_delay_2s/signals.npz")


def test_published_ellipsoid_added_mass_uses_delayed_acceleration():
    with np.load(FIXTURE, allow_pickle=False) as source:
        time = source["time"]
        acceleration = source["delay_input_heave"]
        delayed = source["delay_output_heave"]
        body_acceleration = source["body_acceleration_heave"]
        added_force = source["applied_added_mass_heave"]
        total_force = source["total_hydrodynamic_force_heave"]
        pto_force = source["pto_force_heave"]
        adjusted_mass = float(source["source_adjusted_mass"])
        applied_coefficient = float(
            source["source_applied_added_mass_coefficient"]
        )

    np.testing.assert_allclose(time, np.arange(41) * .05,
                               rtol=0, atol=1e-13)
    np.testing.assert_array_equal(acceleration, body_acceleration)
    np.testing.assert_array_equal(delayed[:2], [0., 0.])
    # The 0.1 microsecond Transport Delay extrapolates from the two prior
    # major-step accelerations. The small delay correction matters when
    # checking this source mechanism at its actual precision.
    expected = (2 * acceleration[1:-1] - acceleration[:-2]
                - (1e-7 / .05) * (acceleration[1:-1]
                                   - acceleration[:-2]))
    assert np.max(np.abs(delayed[2:] - expected)) < 2e-16
    assert np.max(np.abs(added_force - applied_coefficient * delayed)) < 1e-8
    assert np.max(np.abs(adjusted_mass * body_acceleration
                         - total_force - pto_force)) < 1e-7


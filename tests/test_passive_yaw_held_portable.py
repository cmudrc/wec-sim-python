"""Keep held-heading event behavior available in ordinary Python CI."""

import numpy as np

from wecsim.hydrodynamics.passiveYaw import (
    HeldPassiveYawExcitation, SampledPassiveYawExcitation,
)


def test_held_heading_commits_only_accepted_states_and_wraps_at_zero():
    directions = np.arange(0.0, 360.0, 10.0)
    times = np.arange(6) * 0.01
    force_grid = np.zeros((len(times), 1, len(directions), 6))
    force_grid[:, 0, :, 5] = directions
    sampled = SampledPassiveYawExcitation(
        times, directions, np.array([10.0]), force_grid, np.zeros(len(times)),
    )
    held = HeldPassiveYawExcitation(sampled, threshold=1.0)

    # These headings exercise initialization between BEM samples, a hold,
    # snapping to a tabulated heading at the inclusive threshold, an update,
    # and the 0/360-degree boundary. Yaw torque is unchanged by frame rotation.
    yaw_degrees = (2.0, 2.5, 0.0, 1.0, 2.0, 11.0)
    expected_headings = (8.0, 8.0, 10.0, 10.0, 8.0, 0.0)
    for step, (time, yaw_degrees_at_time, expected) in enumerate(zip(
        times, yaw_degrees, expected_headings,
    )):
        coordinate = np.array([np.deg2rad(yaw_degrees_at_time)])
        before = held.last_heading
        trial_force = held(time, coordinate, np.zeros(1))
        assert trial_force[5] == expected
        assert held.last_heading == before
        assert len(held.force_history) == step

        held.commit(time, coordinate)
        assert held.last_heading == expected
        np.testing.assert_allclose(held.force_history[-1], trial_force,
                                   rtol=0, atol=1e-12)

    np.testing.assert_allclose(
        np.asarray(held.force_history)[:, 5], expected_headings,
        rtol=0, atol=1e-12,
    )

"""Analytic cable geometry and tension checks without source traces."""

import numpy as np
import pytest

from wecsim.dynamics.cable import PlanarCableAttachment, WecSimCableTension


def test_slack_cable_and_stretched_cable_follow_signed_source_law():
    cable = WecSimCableTension(stiffness=100, damping=4, length=2,
                              initial_length=1)
    np.testing.assert_array_equal(
        cable.force_z(np.array([0, 1, 1.5]), np.array([2, 2, -3])),
        [0, 0, -38],
    )
    assert cable.force_z(1.5, -20) == 30  # damping is not clamped
    with pytest.raises(ValueError, match="nonnegative"):
        WecSimCableTension(-1, 0, 2, 1)


def test_cable_attachment_rotation_changes_length_and_end_speed():
    cable = PlanarCableAttachment((1, 0), (0, 0), initial_length=1)
    base = np.array([0., 0., 0.])
    follower = np.array([2., 0., 0.])
    still = np.zeros(3)
    assert cable.motion(base, still, follower, still) == (0., 0.)
    stretch, rate = cable.motion(base, [0, 0, 1], follower, still)
    assert stretch == 0
    assert rate == 0  # endpoint initially moves transverse to the cable
    stretch, rate = cable.motion(base, still, [2.5, 0, 0], [2, 0, 0])
    assert stretch == .5
    assert rate == 2
    with pytest.raises(ValueError, match="coincide"):
        cable.motion(base, still, [1, 0, 0], still)


def test_cable_wrenches_conserve_net_force_and_include_attachment_moment():
    cable = PlanarCableAttachment((0, 1), (0, 0), initial_length=1)
    base, follower = [0., 0., 0.], [1., 1., 0.]
    base_wrench, follower_wrench = cable.world_wrenches(base, follower, -20)
    np.testing.assert_allclose(base_wrench, [20, 0, 20])
    np.testing.assert_allclose(follower_wrench, [-20, 0, 0])
    np.testing.assert_allclose(base_wrench[:2] + follower_wrench[:2], 0)
    batch_base = np.tile(base, (2, 1))
    batch_follower = np.array([follower, [0, 2, 0]])
    left, right = cable.world_wrenches(batch_base, batch_follower,
                                       np.array([-20., -10.]))
    assert left.shape == right.shape == (2, 3)
    np.testing.assert_allclose(left[:, :2] + right[:, :2], 0)

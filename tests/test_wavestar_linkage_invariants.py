"""The WaveStar linkage's travel and virtual-work relationships."""

import numpy as np

from wecsim.dynamics.wavestar import WaveStarLinkage


def test_wavestar_arm_motion_is_rotation_about_pivot():
    linkage = WaveStarLinkage()
    angles = np.array([0., .1, -.2])
    positions = linkage.point_position(linkage.pto_arm, angles)
    radius = np.linalg.norm(np.asarray(linkage.pto_arm) - linkage.pivot)
    np.testing.assert_allclose(np.linalg.norm(positions - linkage.pivot, axis=1),
                               radius, rtol=0, atol=1e-15)
    np.testing.assert_allclose(linkage.pto_stroke(0), 0, atol=1e-15)
    speed = linkage.pto_speed(angles, np.ones(3))
    numerical = (linkage.pto_stroke(angles + 1e-6)
                 - linkage.pto_stroke(angles - 1e-6)) / (2e-6)
    np.testing.assert_allclose(speed, numerical, rtol=0, atol=1e-10)
    np.testing.assert_allclose(
        np.sum(linkage.point_velocity(linkage.pto_arm, angles, np.ones(3))
               * (positions - linkage.pivot), axis=1),
        0, atol=1e-15,
    )

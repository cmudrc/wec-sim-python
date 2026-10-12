"""Platform hydro forces on small synthetic histories with known answers."""

import numpy as np
import pytest

from wecsim.hydrodynamics.mostPlatform import MostPlatformHydrodynamics


def _platform():
    return MostPlatformHydrodynamics(
        mass=1e8, inertia=np.diag([1e9, 1e9, 1e9]),
        equilibrium_pose=np.zeros(6),
        restoring_stiffness=np.diag([2., 3., 4., 5., 6., 7.]),
        static_restoring=np.array([1., 2., 3., 4., 5., 6.]),
        drag=np.diag([10., 20., 30., 40., 50., 60.]),
        radiation_time=np.array([0., .1, .2]),
        radiation_kernel=np.repeat(np.eye(6)[None], 3, axis=0),
        added_mass=np.eye(6),
    )


def test_platform_restoring_drag_and_radiation_obey_independent_force_laws():
    platform = _platform()
    pose = np.array([1., -2., 0., 0., .5, 0.])
    np.testing.assert_allclose(platform.restoring_force(pose),
                               [3, -4, 3, 4, 8, 6])
    speed = np.array([2., -3., 0., 0., .5, 0.])
    np.testing.assert_allclose(platform.drag_force(speed),
                               [40, -180, 0, 0, 12.5, 0])
    history = np.tile(speed, (3, 1))
    radiation = platform.radiation_force(history, dt=.1, memory_time=.2)
    np.testing.assert_allclose(radiation[0], 0)
    np.testing.assert_allclose(radiation[1], .1 * speed, atol=1e-12)
    np.testing.assert_allclose(radiation[2], .2 * speed, atol=1e-12)
    with pytest.raises(ValueError, match="valid time steps"):
        platform.radiation_force(history, dt=0)


def test_short_platform_run_responds_to_external_surge_force():
    platform = _platform()
    time = np.arange(4) * .01
    waves = np.zeros((4, 6))
    tower = np.zeros_like(waves)
    reference = platform.simulate(time, waves, tower, radiation_memory=.03)
    waves[:, 0] = 1e7
    driven = platform.simulate(time, waves, tower, radiation_memory=.03)
    assert driven.position.shape == driven.velocity.shape == (4, 6)
    assert np.isfinite(driven.position).all()
    np.testing.assert_allclose(driven.position[0], 0)
    assert driven.position[-1, 0] > reference.position[-1, 0]

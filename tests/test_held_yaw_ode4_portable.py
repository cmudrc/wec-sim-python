"""An analytic motion check for the opt-in held-yaw source solver."""

import numpy as np

from wecsim.dynamics.generalDynamics import (
    BodyMotion, DynamicBody, GeneralizedDynamics,
)
from wecsim.hydrodynamics.passiveYaw import (
    HeldPassiveYawExcitation, SampledPassiveYawExcitation,
)


def test_held_yaw_ode4_constant_torque_and_delayed_heading_commit():
    dt = .1
    directions = np.array([-180., -90., 0., 90., 180.])
    force_grid = np.zeros((5, 1, len(directions), 6))
    force_grid[..., 5] = 2.
    sampled = SampledPassiveYawExcitation(
        np.arange(5) * dt / 2, directions, np.array([0.]),
        force_grid, np.zeros(5),
    )
    held = HeldPassiveYawExcitation(sampled, threshold=1.)
    jacobian = np.zeros((6, 1))
    jacobian[5, 0] = 1.

    def motion(coordinate, speed):
        position = np.zeros(6)
        position[5] = coordinate[0]
        return BodyMotion(position, jacobian, np.zeros(6))

    inertia = np.zeros((6, 6))
    inertia[5, 5] = 1.
    body = DynamicBody(
        rigid_mass=inertia,
        added_mass=(np.zeros((6, 6)),),
        damping=(np.zeros((6, 6)),),
        restoring=np.zeros((6, 6)),
        static_force=np.zeros(6),
        reference_position=np.zeros(6),
        motion=motion,
        excitation=lambda time: np.zeros(6),
        radiation_kernel=np.zeros((3, 6, 6)),
        state_excitation=held,
    )

    result = GeneralizedDynamics((body,), 1).integrate(
        dt=dt, end_time=.2, solver="wecsim_ode4",
    )

    # Unit yaw inertia under 2 Nm has q=t² and v=2t exactly. The model
    # samples excitation at both major and half steps while integrating.
    np.testing.assert_allclose(result.coordinate[:, 0], [0, .01, .04],
                               rtol=0, atol=1e-14)
    np.testing.assert_allclose(result.speed[:, 0], [0, .2, .4],
                               rtol=0, atol=1e-14)
    np.testing.assert_allclose(result.acceleration[:, 0], 2.,
                               rtol=0, atol=1e-14)
    np.testing.assert_allclose(np.asarray(held.force_history)[:, 5], 2.,
                               rtol=0, atol=1e-14)
    assert np.isclose(held.last_heading, 0., rtol=0, atol=1e-12)
    assert held._heading(result.coordinate[-1, 0]) < -2.

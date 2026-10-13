"""Portable numerical guard for the published selected-heading ode4 route."""

import numpy as np

from wecsim.dynamics.generalDynamics import (
    BodyMotion, DynamicBody, GeneralizedDynamics,
)
from wecsim.hydrodynamics.passiveYaw import (
    NearestSampledHeadingExcitation, SampledPassiveYawExcitation,
)


def test_source_ode4_selects_bank_at_rk4_minor_stages():
    """A first-stage heading tie must not hold its force for the whole step."""
    directions = np.arange(0, 360, 10)
    forces = np.zeros((3, 1, len(directions), 6))
    forces[:, 0, 35, 5] = 2  # -10 degree selection
    forces[:, 0, 1, 5] = 4   # +10 degree selection
    sampled = SampledPassiveYawExcitation(
        np.array([0, .05, .1]), directions, np.array([0]),
        forces, np.zeros(3),
    )
    bank = NearestSampledHeadingExcitation(sampled, np.array([-10, 10]))
    jacobian = np.zeros((6, 1))
    jacobian[5, 0] = 1

    def motion(coordinate, speed):
        position = np.zeros(6)
        position[5] = coordinate[0]
        return BodyMotion(position, jacobian, np.zeros(6))

    inertia = np.zeros((6, 6))
    inertia[5, 5] = 1
    body = DynamicBody(
        rigid_mass=inertia,
        added_mass=(np.zeros((6, 6)),),
        damping=(np.zeros((6, 6)),),
        restoring=np.zeros((6, 6)),
        static_force=np.zeros(6),
        reference_position=np.zeros(6),
        motion=motion,
        excitation=lambda time: np.zeros(6),
        radiation_kernel=np.zeros((2, 6, 6)),
        state_excitation=bank,
    )
    result = GeneralizedDynamics((body,), 1).integrate(
        dt=.1, end_time=.1, initial_speed=np.array([-1.0]),
        solver="wecsim_ode4",
    )
    # At q=0, the tie chooses -10 degrees (2 Nm). Negative yaw at the
    # remaining three stages chooses +10 degrees (4 Nm). A major-step hold
    # would end at -0.8 rad/s; stage selection ends at -19/30 rad/s.
    np.testing.assert_allclose(result.coordinate[:, 0],
                               [0, -1 / 12], rtol=0, atol=1e-14)
    np.testing.assert_allclose(result.speed[:, 0],
                               [-1, -19 / 30], rtol=0, atol=1e-14)
    np.testing.assert_allclose(result.acceleration[:, 0],
                               [2, 4], rtol=0, atol=1e-14)



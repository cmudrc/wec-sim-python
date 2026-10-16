"""Independent kinematic and energy checks for the OSWEC pitch hinge."""

import numpy as np
import pytest
from scipy.spatial.transform import Rotation

from wecsim.dynamics.generalDynamics import DynamicBody, GeneralizedDynamics
from wecsim.dynamics.hingePitch import hinged_pitch_motion


@pytest.mark.parametrize("angle,speed", ((0.0, 0.4), (0.37, -0.6), (-0.41, 0.3)))
def test_hinge_jacobian_and_curvature_follow_rigid_rotation(angle, speed):
    lever = 1.7
    motion = hinged_pitch_motion(lever)(np.array([angle]), np.array([speed]))

    def body_displacement(at_time):
        rotated = Rotation.from_euler("y", angle + speed * at_time).apply(
            [0, 0, lever]
        )
        return np.array([
            rotated[0], 0, rotated[2] - lever,
            0, angle + speed * at_time, 0,
        ])

    step = 1e-4
    center = body_displacement(0)
    plus = body_displacement(step)
    minus = body_displacement(-step)
    np.testing.assert_allclose(motion.displacement, center, rtol=0, atol=1e-14)
    np.testing.assert_allclose(motion.jacobian[:, 0] * speed,
                               (plus - minus) / (2 * step), rtol=0, atol=1e-8)
    np.testing.assert_allclose(motion.bias_acceleration,
                               (plus - 2 * center + minus) / step**2,
                               rtol=0, atol=1e-7)


def test_unforced_hinge_loses_only_pto_damper_work():
    """Curved center motion and PTO force must obey one energy balance."""
    mass = np.diag([3.0, 3.0, 3.0, 0.0, 2.0, 0.0])
    added = np.diag([0.4, 0.0, 0.5, 0.0, 0.1, 0.0])
    stiffness, damping = 2.0, 0.8
    body = DynamicBody(
        rigid_mass=mass, added_mass=(added,), damping=(np.zeros((6, 6)),),
        restoring=np.zeros((6, 6)), static_force=np.zeros(6),
        reference_position=np.zeros(6), motion=hinged_pitch_motion(1.7),
        excitation=lambda time: np.zeros(6),
    )
    response = GeneralizedDynamics(
        (body,), 1, pto_stiffness=np.array([[stiffness]]),
        pto_damping=np.array([[damping]]),
    ).integrate(dt=0.005, end_time=5,
                initial_coordinate=np.array([0.4]),
                initial_speed=np.array([0.3]))

    body_speed = response.body_velocity[:, 0]
    kinetic = 0.5 * np.einsum("ti,ij,tj->t", body_speed, mass + added, body_speed)
    spring = 0.5 * stiffness * response.coordinate[:, 0]**2
    dissipated = damping * np.r_[
        0, np.cumsum(np.diff(response.time)
                     * (response.speed[:-1, 0]**2
                        + response.speed[1:, 0]**2) / 2),
    ]
    total = kinetic + spring + dissipated
    assert np.ptp(response.coordinate[:, 0]) > 0.2
    assert dissipated[-1] > 0.1
    assert np.max(np.abs(total - total[0])) / total[0] < 2e-5

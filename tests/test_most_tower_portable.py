"""Check tower reaction mechanics without the large MATLAB MOST artifacts."""

import numpy as np

from wecsim import MostTowerReaction


def test_collocated_turbine_reaction_obeys_newton_euler_balance():
    # Put every mass at the tower base. This removes lever-arm moments while
    # retaining the tower, yaw bearing, nacelle, hub, and three blade paths.
    origin = np.zeros(3)
    properties = {
        "tower": {"offset": 0, "mass": 100, "cog": origin,
                  "Inertia": origin},
        "nacelle": {"mass_yawBearing": 2, "cog_yawBearing": origin,
                    "mass": 3, "cog": origin, "Inertia": origin,
                    "tiltangle": 0},
        "hub": {"mass": 4, "cog": origin, "Rhub": 0, "precone": 0,
                "Inertia": origin},
        "blade": {"mass": 5, "cog_rel": origin, "Inertia": origin},
        "Inertia_Rotor_cogHub": 1,
    }
    tower = MostTowerReaction(properties, origin)
    total_mass = 100 + 2 + 3 + 4 + 3 * 5
    position = np.zeros((2, 6))
    velocity = np.zeros_like(position)
    acceleration = np.zeros_like(position)
    acceleration[1, 2] = 2
    rotor_speed = np.zeros(2)
    azimuth = np.array([0, .7])
    torque = np.zeros(2)
    blade_loads = np.zeros((2, 6, 3))

    reaction = tower.evaluate(position, velocity, acceleration,
                              rotor_speed, azimuth, torque, blade_loads)
    np.testing.assert_allclose(reaction[:, :2], 0, atol=1e-12)
    np.testing.assert_allclose(reaction[:, 3:], 0, atol=1e-12)
    np.testing.assert_allclose(
        reaction[:, 2], -total_mass * (tower.gravity + np.array([0, 2])),
        rtol=0, atol=1e-12,
    )

    # An axial force at each blade root adds three newtons at the base,
    # regardless of azimuth. The mass reaction and gravity remain separate.
    blade_loads[:, 0, :] = 1
    driven = tower.evaluate(position, velocity, acceleration,
                            rotor_speed, azimuth, torque, blade_loads)
    np.testing.assert_allclose(driven[:, 0] - reaction[:, 0], 3,
                               rtol=0, atol=1e-12)
    np.testing.assert_allclose(driven[:, 1:] - reaction[:, 1:], 0,
                               rtol=0, atol=1e-12)

    matrix = tower.acceleration_matrix(azimuth)
    np.testing.assert_allclose(matrix[:, :3, :3],
                               np.tile(total_mass * np.eye(3), (2, 1, 1)),
                               rtol=0, atol=1e-12)

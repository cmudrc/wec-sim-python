"""Tower-base reaction of the published MOST IEA 15 MW turbine."""

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.io import loadmat

from wecsim.reference_models.mostBEM import _cross3, _rx, _ry
from wecsim.reference_models.mostMooring import MostStaticMooring


@dataclass(frozen=True)
class MostTowerReaction:
    """Newton–Euler reaction on the platform from the MOST turbine.

    Coordinates and velocities are world-frame platform-center records in
    surge/sway/heave/roll/pitch/yaw order. Root loads use the three preconed
    blade frames returned by :class:`MostBEM`. The returned wrench is in the
    rotating platform frame, about the tower base, as in WEC-Sim's MOST log.
    """

    properties: dict
    platform_cg: np.ndarray
    gravity: float = 9.80665

    @classmethod
    def from_iea15mw(cls, properties_file: str | Path, *,
                     platform_cg=(0, 0, -14.4)) -> "MostTowerReaction":
        """Read the generated ``Properties_IEA15MW.mat`` source artifact."""
        properties = loadmat(properties_file, simplify_cells=True)["WTcomponents"]
        cg = np.asarray(platform_cg, dtype=float)
        if cg.shape != (3,) or not np.isfinite(cg).all():
            raise ValueError("platform center of gravity must be a finite three-vector")
        return cls(properties, cg)

    def acceleration_matrix(self, azimuth) -> np.ndarray:
        """Return N×6×6 turbine inertia about the tower base.

        The matrix maps world-frame platform-center acceleration to the
        platform-frame tower reaction: ``load = zero_acceleration_load - M @ a``.
        It depends on blade azimuth, but not wind, rotor speed, or torque.
        """
        azimuth = np.asarray(azimuth, dtype=float).reshape(-1)
        if azimuth.size < 1 or not np.isfinite(azimuth).all():
            raise ValueError("azimuth must contain finite values")
        n = azimuth.size
        position = np.tile(np.r_[self.platform_cg, np.zeros(3)], (n, 1))
        zero_state = np.zeros((n, 6))
        zero_scalar = np.zeros(n)
        zero_load = np.zeros((n, 6, 3))
        bias = self.evaluate(position, zero_state, zero_state, zero_scalar,
                             azimuth, zero_scalar, zero_load)
        matrix = np.empty((n, 6, 6))
        for axis in range(6):
            acceleration = np.zeros((n, 6))
            acceleration[:, axis] = 1
            matrix[:, :, axis] = bias - self.evaluate(
                position, zero_state, acceleration, zero_scalar, azimuth,
                zero_scalar, zero_load,
            )
        return matrix

    def evaluate(self, position, velocity, acceleration, rotor_speed, azimuth,
                 generator_torque, blade_root_load) -> np.ndarray:
        """Return N×6 tower-base loads for aligned turbine/platform states.

        ``rotor_speed`` is in rad/s, ``azimuth`` in rad, generator torque in
        N m, and ``blade_root_load`` has shape N×6×3. Platform acceleration
        includes both linear and angular components. Rotor acceleration is
        obtained from the blade-root shaft torque and generator torque.
        """
        position = np.asarray(position, dtype=float)
        velocity = np.asarray(velocity, dtype=float)
        acceleration = np.asarray(acceleration, dtype=float)
        n = len(position) if position.ndim == 2 else 0
        speed = np.asarray(rotor_speed, dtype=float).reshape(-1)
        azimuth = np.asarray(azimuth, dtype=float).reshape(-1)
        torque = np.asarray(generator_torque, dtype=float).reshape(-1)
        loads = np.asarray(blade_root_load, dtype=float)
        if (n < 1 or any(x.shape != (n, 6) for x in
                         (position, velocity, acceleration))
                or any(x.shape != (n,) for x in (speed, azimuth, torque))
                or loads.shape != (n, 6, 3)
                or not all(np.isfinite(x).all() for x in
                           (position, velocity, acceleration, speed,
                            azimuth, torque, loads))):
            raise ValueError("MOST tower needs aligned finite six-DOF states and three blade loads")

        p = self.properties
        tower, nacelle, hub, blade = (p[name] for name in
                                      ("tower", "nacelle", "hub", "blade"))
        base = np.array([0., 0., tower["offset"]])
        hub_pos = np.asarray(hub["cog"], dtype=float)
        hub_radius = float(hub["Rhub"])
        cone = np.deg2rad(hub["precone"])
        tilt = _ry(np.deg2rad(nacelle["tiltangle"]))
        precone = _ry(-cone)
        blade_arm = np.array([0., 0., hub_radius + blade["cog_rel"][2]])
        root_arm = np.array([0., 0., hub_radius])
        gravity = np.array([0., 0., -self.gravity])
        fixed = (
            (tower["mass"], tower["cog"], np.diag(tower["Inertia"])),
            (nacelle["mass_yawBearing"], nacelle["cog_yawBearing"],
             np.zeros((3, 3))),
            (nacelle["mass"], nacelle["cog"], np.diag(nacelle["Inertia"])),
        )
        hub_inertia = np.diag(hub["Inertia"])
        blade_inertia = np.diag(blade["Inertia"])
        result = np.empty((n, 6))

        for i in range(n):
            rotation = MostStaticMooring._rotation(*position[i, 3:])
            omega = velocity[i, 3:]
            alpha = acceleration[i, 3:]
            cg_acc = acceleration[i, :3]
            inertial_force = np.zeros(3)
            inertial_moment = np.zeros(3)
            gravity_force = np.zeros(3)
            gravity_moment = np.zeros(3)

            def add_component(mass, local_pos, inertia, angular_speed,
                              angular_acceleration, linear_acceleration):
                nonlocal inertial_force, inertial_moment
                nonlocal gravity_force, gravity_moment
                arm = rotation @ (np.asarray(local_pos) - base)
                force = mass * linear_acceleration
                weight = mass * gravity
                inertial_force += force
                inertial_moment += (_cross3(arm, force)
                                    + inertia @ angular_acceleration
                                    + _cross3(angular_speed,
                                               inertia @ angular_speed))
                gravity_force += weight
                gravity_moment += _cross3(arm, weight)

            def fixed_acceleration(local_pos):
                arm = rotation @ (np.asarray(local_pos) - self.platform_cg)
                return (cg_acc + _cross3(alpha, arm)
                        + _cross3(omega, _cross3(omega, arm)))

            for mass, local_pos, inertia in fixed:
                add_component(mass, local_pos, rotation @ inertia @ rotation.T,
                              omega, alpha, fixed_acceleration(local_pos))

            shaft = rotation @ tilt @ np.array([1., 0., 0.])
            spin = speed[i] * shaft
            aero_shaft = np.sum(
                np.cos(cone)*loads[i, 3] - np.sin(cone)*loads[i, 5]
                - hub_radius*np.cos(cone)*loads[i, 1]
            )
            rotor_acc = ((aero_shaft - torque[i]) / p["Inertia_Rotor_cogHub"]
                         - np.dot(alpha, shaft))
            rotor_omega = omega + spin
            rotor_alpha = alpha + rotor_acc*shaft + _cross3(omega, spin)
            hub_acc = fixed_acceleration(hub_pos)
            hub_rotation = rotation @ tilt
            add_component(hub["mass"], hub_pos,
                          hub_rotation @ hub_inertia @ hub_rotation.T,
                          rotor_omega, rotor_alpha, hub_acc)

            aero_force = np.zeros(3)
            aero_moment = np.zeros(3)
            for blade_index in range(3):
                local_root = tilt @ _rx(azimuth[i] + 2*np.pi*blade_index/3) @ precone
                world_root = rotation @ local_root
                relative_blade_pos = world_root @ blade_arm
                blade_pos = hub_pos + local_root @ blade_arm
                blade_acc = (hub_acc
                             + _cross3(rotor_alpha, relative_blade_pos)
                             + _cross3(rotor_omega,
                                       _cross3(rotor_omega, relative_blade_pos)))
                add_component(blade["mass"], blade_pos,
                              world_root @ blade_inertia @ world_root.T,
                              rotor_omega, rotor_alpha, blade_acc)
                force = world_root @ loads[i, :3, blade_index]
                moment = world_root @ loads[i, 3:, blade_index]
                root_pos = hub_pos + local_root @ root_arm
                aero_force += force
                aero_moment += moment + _cross3(rotation @ (root_pos - base),
                                                force)

            result[i, :3] = rotation.T @ (
                aero_force + gravity_force - inertial_force
            )
            result[i, 3:] = rotation.T @ (
                aero_moment + gravity_moment - inertial_moment
            )
        return result

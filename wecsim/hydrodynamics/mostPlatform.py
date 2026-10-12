"""Hydrodynamic force laws for the published MOST VolturnUS platform."""

from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np
from scipy.interpolate import CubicSpline

from wecsim.dynamics.generalDynamics import BodyMotion, DynamicBody, GeneralizedDynamics
from wecsim.dynamics.mostMooring import MostStaticMooring
from wecsim.dynamics.mostTower import MostTowerReaction


@dataclass(frozen=True)
class MostPlatformResponse:
    time: np.ndarray
    position: np.ndarray  # world-frame center-of-gravity pose, m and rad
    velocity: np.ndarray  # world-frame center-of-gravity velocity, m/s and rad/s
    acceleration: np.ndarray  # active independent-coordinate acceleration


@dataclass(frozen=True)
class MostPlatformHydrodynamics:
    """Evaluate source-case restoring, drag, and radiation on platform motion.

    Forces follow WEC-Sim's logged resisting-force convention. ``pose`` is
    the six-component position of the platform center of gravity in world
    coordinates. ``simulate`` advances the dominant surge, heave, and pitch
    coordinates when a tower-base reaction history is supplied.
    """

    mass: float
    inertia: np.ndarray
    equilibrium_pose: np.ndarray
    restoring_stiffness: np.ndarray
    static_restoring: np.ndarray
    drag: np.ndarray
    radiation_time: np.ndarray
    radiation_kernel: np.ndarray
    added_mass: np.ndarray
    tower_base_height: float = 15.0  # world height at zero platform motion

    @classmethod
    def from_volturnus(cls, h5_file: str | Path,
                       mass_properties_file: str | Path) -> "MostPlatformHydrodynamics":
        """Read the pinned MOST hydrodynamics and platform mass properties."""
        rho, gravity = 1025.0, 9.80665
        with h5py.File(h5_file) as hydro:
            body = hydro["body1"]
            volume = float(body["properties/disp_vol"][0, 0])
            center_gravity = body["properties/cg"][:].ravel()
            center_buoyancy = body["properties/cb"][:].ravel()
            coefficients = body["hydro_coeffs"]
            stiffness = coefficients["linear_restoring_stiffness"][:] * rho * gravity
            added_mass = coefficients["added_mass/inf_freq"][:] * rho
            irf = coefficients["radiation_damping/impulse_response_fun"]
            radiation_time = irf["t"][:].ravel()
            radiation_kernel = np.moveaxis(irf["K"][:], -1, 0) * rho
        with h5py.File(mass_properties_file) as properties:
            platform = properties["Platform/VolturnUS"]
            mass = float(platform["mass"][0, 0])
            inertia = platform["I_COG"][:]

        buoyancy = rho * gravity * volume
        restoring = np.zeros(6)
        restoring[2] = mass * gravity - buoyancy
        restoring[3:] = -np.cross(
            center_buoyancy - center_gravity, [0, 0, buoyancy],
        )
        drag = np.array([
            [9.23e5, 0, 0, 0, -8.92e6, 0],
            [0, 9.23e5, 0, 8.92e6, 0, 0],
            [0, 0, 2.30e6, 0, 0, 0],
            [0, 8.92e6, 0, 1.68e10, 0, 0],
            [-8.92e6, 0, 0, 0, 1.68e10, 0],
            [0, 0, 0, 0, 0, 4.80e10],
        ])
        return cls(
            mass=mass, inertia=inertia,
            equilibrium_pose=np.r_[center_gravity, np.zeros(3)],
            restoring_stiffness=stiffness, static_restoring=restoring,
            drag=drag, radiation_time=radiation_time,
            radiation_kernel=radiation_kernel, added_mass=added_mass,
        )

    def restoring_force(self, pose: np.ndarray) -> np.ndarray:
        """Return the logged resisting hydrostatic force on saved or live poses."""
        pose = np.asarray(pose, dtype=float)
        if pose.shape[-1:] != (6,) or not np.isfinite(pose).all():
            raise ValueError("platform pose must have six finite coordinates")
        return ((pose - self.equilibrium_pose) @ self.restoring_stiffness.T
                + self.static_restoring)

    def drag_force(self, velocity: np.ndarray) -> np.ndarray:
        """Return the logged quadratic viscous force on saved or live speeds."""
        velocity = np.asarray(velocity, dtype=float)
        if velocity.shape[-1:] != (6,) or not np.isfinite(velocity).all():
            raise ValueError("platform velocity must have six finite coordinates")
        return (velocity * np.abs(velocity)) @ self.drag.T

    def radiation_force(self, velocity: np.ndarray, dt: float,
                        memory_time: float = 60.0) -> np.ndarray:
        """Convolve platform velocity with the WEC-Sim radiation kernel."""
        velocity = np.asarray(velocity, dtype=float)
        if (velocity.ndim != 2 or velocity.shape[1] != 6
                or not np.isfinite(velocity).all()
                or not np.isfinite([dt, memory_time]).all()
                or dt <= 0 or memory_time <= 0
                or memory_time > self.radiation_time[-1]):
            raise ValueError("radiation needs finite six-DOF history and valid time steps")
        lags = min(len(velocity) - 1, int(round(memory_time / dt)))
        kernel = CubicSpline(self.radiation_time, self.radiation_kernel,
                             axis=0)(np.arange(lags + 1) * dt)
        force = np.zeros_like(velocity)
        for index in range(1, len(velocity)):
            last = min(index, lags)
            force[index] = dt * np.einsum(
                "tij,tj->i", kernel[:last + 1],
                velocity[index - last:index + 1][::-1],
            )
            force[index] -= dt / 2 * (
                kernel[0] @ velocity[index]
                + kernel[last] @ velocity[index - last]
            )
        return force

    def simulate(self, time: np.ndarray, wave_excitation: np.ndarray,
                 tower_base_load: np.ndarray, *,
                 mooring: MostStaticMooring | None = None,
                 radiation_memory: float = 60.0) -> MostPlatformResponse:
        """Advance published surge, heave, and pitch with a prescribed tower load.

        ``wave_excitation`` is a Python-generated six-component HDF5 force.
        ``tower_base_load`` is a six-component force/moment history in the
        platform frame. The tower history is an external input; this method
        does not model turbine-to-platform feedback or the three inactive
        platform coordinates.
        """
        return self._simulate(time, wave_excitation, tower_base_load,
                              mooring=mooring, radiation_memory=radiation_memory)

    def simulate_with_turbine(
            self, time: np.ndarray, wave_excitation: np.ndarray,
            tower: MostTowerReaction, rotor_speed: np.ndarray,
            azimuth: np.ndarray, generator_torque: np.ndarray,
            blade_root_load: np.ndarray, *,
            mooring: MostStaticMooring | None = None,
            radiation_memory: float = 60.0,
            full_six_dof: bool = True,
            turbine_step=None) -> MostPlatformResponse:
        """Advance platform motion with the turbine's live inertia and force.

        Rotor state, generator torque, and blade-root loads are prescribed
        histories; the tower reaction is recomputed at each platform state.
        This does not independently advance the rotor. All six platform
        coordinates are active by default; the reduced three-coordinate
        comparison remains available with ``full_six_dof=False``. A
        ``turbine_step`` callback advances mutable rotor histories once from
        the committed previous state and predicted next platform state.
        """
        if not isinstance(tower, MostTowerReaction):
            raise TypeError("tower must be MostTowerReaction")
        if not np.allclose(tower.platform_cg, self.equilibrium_pose[:3],
                           rtol=0, atol=1e-9):
            raise ValueError("tower and platform centers of gravity differ")
        if not np.isclose(tower.properties["tower"]["offset"],
                          self.tower_base_height, rtol=0, atol=1e-9):
            raise ValueError("tower-base heights differ")
        time = np.asarray(time, dtype=float)
        n = len(time) if time.ndim == 1 else 0
        rotor_speed = np.asarray(rotor_speed, dtype=float).reshape(-1)
        azimuth = np.asarray(azimuth, dtype=float).reshape(-1)
        generator_torque = np.asarray(generator_torque, dtype=float).reshape(-1)
        blade_root_load = np.asarray(blade_root_load, dtype=float)
        if (n < 2 or any(x.shape != (n,) for x in
                         (rotor_speed, azimuth, generator_torque))
                or blade_root_load.shape != (n, 6, 3)
                or not all(np.isfinite(x).all() for x in
                           (rotor_speed, azimuth, generator_torque,
                            blade_root_load))):
            raise ValueError("MOST turbine needs aligned finite rotor and blade histories")
        if turbine_step is not None and not callable(turbine_step):
            raise TypeError("turbine_step must be callable")
        return self._simulate(
            time, wave_excitation, None, mooring=mooring,
            radiation_memory=radiation_memory,
            turbine=(tower, rotor_speed, azimuth, generator_torque,
                     blade_root_load),
            coordinate_axes=tuple(range(6)) if full_six_dof else (0, 2, 4),
            turbine_step=turbine_step,
        )

    def _simulate(self, time, wave_excitation, tower_base_load, *, mooring,
                  radiation_memory, turbine=None,
                  coordinate_axes=(0, 2, 4), turbine_step=None) -> MostPlatformResponse:
        time = np.asarray(time, dtype=float)
        wave_excitation = np.asarray(wave_excitation, dtype=float)
        if turbine is None:
            tower_base_load = np.asarray(tower_base_load, dtype=float)
        if (time.ndim != 1 or len(time) < 2 or not np.isfinite(time).all()
                or not np.isclose(time[0], 0, rtol=0, atol=1e-12)
                or wave_excitation.shape != (len(time), 6)
                or not np.isfinite(wave_excitation).all()
                or (turbine is None and (
                    tower_base_load.shape != (len(time), 6)
                    or not np.isfinite(tower_base_load).all()))):
            raise ValueError("MOST platform needs aligned finite time and six-component loads")
        dt = time[1] - time[0]
        if (dt <= 0 or not np.allclose(np.diff(time), dt, rtol=0, atol=1e-12)
                or not np.isfinite(radiation_memory)
                or radiation_memory <= 0):
            raise ValueError("MOST platform needs uniform positive time and radiation memory")
        memory = min(radiation_memory, time[-1])
        if memory > self.radiation_time[-1]:
            raise ValueError("radiation memory exceeds the hydrodynamic IRF")
        kernel = CubicSpline(self.radiation_time, self.radiation_kernel,
                             axis=0)(np.arange(round(memory / dt) + 1) * dt)
        mooring = MostStaticMooring() if mooring is None else mooring
        if not isinstance(mooring, MostStaticMooring):
            raise TypeError("mooring must be MostStaticMooring")

        projection = np.eye(6)[:, coordinate_axes]
        mooring_offset = -self.equilibrium_pose[:3]
        tower_offset = mooring_offset + np.array([0, 0, self.tower_base_height])
        if turbine is not None:
            tower, rotor_speed, azimuth, generator_torque, blade_root_load = turbine
            if turbine_step is None:
                local_mass = tower.acceleration_matrix(azimuth)
            else:
                local_mass = np.zeros((len(time), 6, 6))
                local_mass[0] = tower.acceleration_matrix(azimuth[:1])[0]
            zero_acceleration = np.zeros((1, 6))

        def motion(coordinate, speed):
            displacement = projection @ coordinate
            rates = projection @ speed
            _, pitch, yaw = displacement[3:]
            # The pose uses roll/pitch/yaw; WEC-Sim reports world angular
            # velocity, which is not the derivative of those three angles.
            cp, sp = np.cos(pitch), np.sin(pitch)
            cy, sy = np.cos(yaw), np.sin(yaw)
            angular_jacobian = np.array([
                [cy * cp, -sy, 0],
                [sy * cp, cy, 0],
                [-sp, 0, 1],
            ])
            jacobian = projection.copy()
            jacobian[3:] = angular_jacobian @ projection[3:]
            roll_rate, pitch_rate, yaw_rate = rates[3:]
            angular_bias = (
                roll_rate * pitch_rate * np.array([-cy * sp, -sy * sp, -cp])
                + roll_rate * yaw_rate * np.array([-sy * cp, cy * cp, 0])
                + pitch_rate * yaw_rate * np.array([-cy, -sy, 0])
            )
            return BodyMotion(displacement, jacobian,
                              np.r_[np.zeros(3), angular_bias])

        def state(at_time, coordinate, speed):
            index = min(len(time) - 1, max(0, int(round(at_time / dt))))
            body_motion = motion(coordinate, speed)
            pose = self.equilibrium_pose + body_motion.displacement
            velocity = body_motion.jacobian @ speed
            rotation = mooring._rotation(*pose[3:])
            arm_mooring = rotation @ mooring_offset
            arm_tower = rotation @ tower_offset
            return index, pose, velocity, rotation, arm_mooring, arm_tower

        def live_force(at_time, coordinate, speed):
            index, pose, velocity, rotation, arm_mooring, arm_tower = state(
                at_time, coordinate, speed,
            )
            mooring_pose = np.r_[pose[:3] + arm_mooring, pose[3:]]
            mooring_load = mooring.force(mooring_pose)[0]
            mooring_wrench = np.r_[
                mooring_load[:3],
                mooring_load[3:] + np.cross(arm_mooring, mooring_load[:3]),
            ]
            force = (wave_excitation[index] + mooring_wrench
                     - self.restoring_force(pose) - self.drag_force(velocity))
            if turbine is None:
                tower_load = tower_base_load[index]
                tower_force = rotation @ tower_load[:3]
                force += np.r_[
                    tower_force,
                    rotation @ tower_load[3:]
                    + np.cross(arm_tower, tower_force),
                ]
            return force

        step_force = None
        if turbine_step is not None:
            previous_coordinate = np.zeros(len(coordinate_axes))
            previous_speed = np.zeros(len(coordinate_axes))

            def commit_state(_at_time, coordinate, speed):
                previous_coordinate[:] = coordinate
                previous_speed[:] = speed

            live_force.commit_state = commit_state

            def body_state(coordinate, speed):
                body_motion = motion(coordinate, speed)
                return (self.equilibrium_pose + body_motion.displacement,
                        body_motion.jacobian @ speed)

            def step_force(at_time, interval, predicted_coordinate,
                           predicted_speed):
                index = int(round((at_time + interval) / dt))
                previous_pose, previous_velocity = body_state(
                    previous_coordinate, previous_speed,
                )
                predicted_pose, predicted_velocity = body_state(
                    predicted_coordinate, predicted_speed,
                )
                turbine_step(
                    index, at_time, interval,
                    previous_pose, previous_velocity,
                    predicted_pose, predicted_velocity,
                )
                local_mass[index] = tower.acceleration_matrix(
                    azimuth[index:index + 1],
                )[0]
                return np.zeros(len(coordinate_axes))

        def tower_inertia(at_time, coordinate, speed):
            index, pose, velocity, rotation, _, arm_tower = state(
                at_time, coordinate, speed,
            )
            bias = tower.evaluate(
                pose[None], velocity[None], zero_acceleration,
                rotor_speed[index:index + 1], azimuth[index:index + 1],
                generator_torque[index:index + 1],
                blade_root_load[index:index + 1],
            )[0]
            arm_cross = np.array([
                [0, -arm_tower[2], arm_tower[1]],
                [arm_tower[2], 0, -arm_tower[0]],
                [-arm_tower[1], arm_tower[0], 0],
            ])
            transform = np.block([
                [rotation, np.zeros((3, 3))],
                [arm_cross @ rotation, rotation],
            ])
            to_local = np.block([
                [rotation.T, np.zeros((3, 3))],
                [np.zeros((3, 3)), rotation.T],
            ])
            # Simscape rotates rigid inertia with the platform and includes
            # its gyroscopic moment; the hydro added mass remains separate.
            inertia_world = rotation @ self.inertia @ rotation.T
            platform_force = np.r_[
                np.zeros(3),
                -np.cross(velocity[3:], inertia_world @ velocity[3:]),
            ]
            platform_mass = np.zeros((6, 6))
            platform_mass[3:, 3:] = inertia_world - self.inertia
            return (transform @ bias + platform_force,
                    transform @ local_mass[index] @ to_local + platform_mass)

        rigid_mass = np.diag([self.mass] * 3 + [0.] * 3)
        rigid_mass[3:, 3:] = self.inertia
        body = DynamicBody(
            rigid_mass=rigid_mass, added_mass=(self.added_mass,),
            damping=(np.zeros((6, 6)),), restoring=np.zeros((6, 6)),
            static_force=np.zeros(6), reference_position=self.equilibrium_pose,
            motion=motion, excitation=lambda _: np.zeros(6),
            radiation_kernel=kernel, state_excitation=live_force,
            state_inertia=tower_inertia if turbine is not None else None,
        )
        response = GeneralizedDynamics(
            (body,), len(coordinate_axes), added_mass_delay=1e-7,
        ).integrate(dt=dt, end_time=float(time[-1]), step_force=step_force)
        return MostPlatformResponse(
            response.time, response.body_position[:, 0],
            response.body_velocity[:, 0], response.acceleration,
        )

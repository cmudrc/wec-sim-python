"""Bidirectionally coupled MOST platform and turbine motion."""

from dataclasses import dataclass

import numpy as np

from wecsim.reference_models.mostMooring import MostStaticMooring
from wecsim.reference_models.mostPlatform import MostPlatformHydrodynamics, MostPlatformResponse
from wecsim.reference_models.mostRotor import MostRotor, MostRotorResponse
from wecsim.reference_models.mostTower import MostTowerReaction
from wecsim.environment.turbSim import MostConstantWind, MostWindField


@dataclass(frozen=True)
class MostCoupledResponse:
    platform: MostPlatformResponse
    rotor: MostRotorResponse
    iterations: int  # one causal pass or the number of global iterations
    position_residual: float  # m or rad, predictor or last-pass correction
    velocity_residual: float  # m/s or rad/s, predictor or last-pass correction


@dataclass(frozen=True)
class MostCoupled:
    """Advance the published MOST turbine with six-coordinate platform motion.

    Each pass advances the Python rotor from the current platform trajectory,
    then solves platform motion with the resulting blade loads and implicit
    turbine inertia. Passes continue until both trajectories are consistent.
    ``full_six_dof=False`` retains the earlier surge/heave/pitch reduction.
    """

    platform: MostPlatformHydrodynamics
    rotor: MostRotor
    tower: MostTowerReaction

    def simulate_causal(self, time, wave_excitation,
                        wind: MostWindField | MostConstantWind,
                        *, mooring: MostStaticMooring | None = None
                        ) -> MostCoupledResponse:
        """Advance turbine and platform once per step using a predicted pose.

        The turbine state is committed before the platform's implicit solve
        for that step. Residuals report the largest predictor correction.
        """
        if (not isinstance(self.platform, MostPlatformHydrodynamics)
                or not isinstance(self.rotor, MostRotor)
                or not isinstance(self.tower, MostTowerReaction)
                or not isinstance(wind, (MostWindField, MostConstantWind))):
            raise TypeError("MOST causal coupling needs platform, rotor, tower, and wind")
        time = np.asarray(time, dtype=float)
        if (time.ndim != 1 or time.size < 2 or not np.isfinite(time).all()
                or not np.all(np.diff(time) > 0)):
            raise ValueError("MOST causal coupling needs increasing finite time")
        count = time.size
        rotor = self.rotor
        speed = np.zeros(count)
        azimuth = np.zeros(count)
        torque = np.zeros(count)
        pitch = np.zeros(count)
        blade_load = np.zeros((count, 6, 3))
        predicted_body_position = np.zeros((count, 6))
        predicted_body_velocity = np.zeros((count, 6))
        speed[0] = rotor.controller.initial_omega
        azimuth[0] = rotor.initial_azimuth
        controller_state = np.array([0., speed[0] / 5, 0.])
        torque[0], pitch[0] = rotor.controller._command(
            controller_state, rotor.controller.initial_pitch,
        )
        precone = np.deg2rad(4)
        cone_cosine, cone_sine = np.cos(precone), np.sin(precone)

        def aerodynamic_load(at_time, pose, velocity, rotor_speed,
                             rotor_azimuth, blade_pitch):
            rotation = MostStaticMooring._rotation(*pose[3:])
            arm = rotation @ rotor.hub_offset
            hub_state = np.r_[
                pose[:3] + arm, pose[3:],
                velocity[:3] + np.cross(velocity[3:], arm),
                velocity[3:], rotor_azimuth, rotor_speed,
            ]
            loads = rotor.bem.loads(
                hub_state, blade_pitch, wind.sampler_at(float(at_time)),
            )
            if not np.isfinite(loads).all():
                raise ValueError("MOST BEM loads are not finite at the platform state")
            shaft_speed = float(velocity[3:] @ (rotation @ rotor.shaft_axis))
            return loads, shaft_speed

        def shaft_torque(loads):
            return float(np.sum(
                cone_cosine * loads[3] - cone_sine * loads[5]
                - rotor.bem.hub_radius * cone_cosine * loads[1]
            ))

        initial_pose = self.platform.equilibrium_pose
        predicted_body_position[0] = initial_pose
        blade_load[0], _ = aerodynamic_load(
            time[0], initial_pose, np.zeros(6),
            speed[0], azimuth[0], pitch[0],
        )

        def turbine_step(index, previous_time, interval,
                         previous_pose, previous_velocity,
                         predicted_pose, predicted_velocity):
            nonlocal controller_state
            predicted_body_position[index] = predicted_pose
            predicted_body_velocity[index] = predicted_velocity
            previous_load, previous_shaft_speed = aerodynamic_load(
                previous_time, previous_pose, previous_velocity,
                speed[index - 1], azimuth[index - 1], pitch[index - 1],
            )
            predicted_shaft_speed = float(
                predicted_velocity[3:] @ (
                    MostStaticMooring._rotation(*predicted_pose[3:])
                    @ rotor.shaft_axis
                )
            )
            # Use the shaft-speed change over this step; the future
            # platform state needed by a centered gradient is unavailable.
            shaft_acceleration = (
                predicted_shaft_speed - previous_shaft_speed
            ) / interval
            first_acceleration = (
                (shaft_torque(previous_load) - torque[index - 1])
                / rotor.inertia - shaft_acceleration
            )
            predicted_speed = speed[index - 1] + interval * first_acceleration
            predicted_azimuth = azimuth[index - 1] + interval * speed[index - 1]
            predicted_load, _ = aerodynamic_load(
                time[index], predicted_pose, predicted_velocity,
                predicted_speed, predicted_azimuth, pitch[index - 1],
            )
            second_acceleration = (
                (shaft_torque(predicted_load) - torque[index - 1])
                / rotor.inertia - shaft_acceleration
            )
            speed[index] = speed[index - 1] + interval * (
                first_acceleration + second_acceleration
            ) / 2
            azimuth[index] = azimuth[index - 1] + interval * (
                speed[index - 1] + speed[index]
            ) / 2
            controller_state, torque[index], pitch[index] = (
                rotor.controller._advance(
                    controller_state, speed[index - 1], speed[index], interval,
                    torque[index - 1], pitch[index - 1],
                )
            )
            blade_load[index], _ = aerodynamic_load(
                time[index], predicted_pose, predicted_velocity,
                speed[index], azimuth[index], pitch[index],
            )

        platform_response = self.platform.simulate_with_turbine(
            time, wave_excitation, self.tower,
            speed, azimuth, torque, blade_load,
            mooring=mooring, full_six_dof=True,
            turbine_step=turbine_step,
        )
        rotor_response = MostRotorResponse(
            time, speed, azimuth, torque, pitch, blade_load,
        )
        return MostCoupledResponse(
            platform_response, rotor_response, 1,
            float(np.max(np.abs(
                platform_response.position - predicted_body_position
            ))),
            float(np.max(np.abs(
                platform_response.velocity - predicted_body_velocity
            ))),
        )

    def simulate(self, time, wave_excitation,
                 wind: MostWindField | MostConstantWind, *,
                 position_tolerance: float = 1e-6,
                 velocity_tolerance: float = 1e-6,
                 max_iterations: int = 20,
                 mooring: MostStaticMooring | None = None,
                 full_six_dof: bool = True) -> MostCoupledResponse:
        """Use only wave and wind inputs to advance the coupled case."""
        if (not isinstance(self.platform, MostPlatformHydrodynamics)
                or not isinstance(self.rotor, MostRotor)
                or not isinstance(self.tower, MostTowerReaction)):
            raise TypeError("MOST coupling needs platform, rotor, and tower models")
        if (not np.isfinite([position_tolerance, velocity_tolerance]).all()
                or position_tolerance <= 0 or velocity_tolerance <= 0
                or not isinstance(max_iterations, int)
                or max_iterations < 1):
            raise ValueError("MOST coupling needs positive tolerances and iterations")
        time = np.asarray(time, dtype=float)
        if time.ndim != 1 or time.size < 2:
            raise ValueError("MOST coupling needs a time history")
        position = np.tile(self.platform.equilibrium_pose, (time.size, 1))
        velocity = np.zeros_like(position)
        axes = tuple(range(6)) if full_six_dof else (0, 2, 4)
        for iteration in range(1, max_iterations + 1):
            rotor_response = self.rotor.simulate(
                time, position, velocity, wind,
            )
            platform_response = self.platform.simulate_with_turbine(
                time, wave_excitation, self.tower,
                rotor_response.rotor_speed, rotor_response.azimuth,
                rotor_response.generator_torque,
                rotor_response.blade_root_load, mooring=mooring,
                full_six_dof=full_six_dof,
            )
            position_residual = float(np.max(np.abs(
                platform_response.position[:, axes] - position[:, axes]
            )))
            velocity_residual = float(np.max(np.abs(
                platform_response.velocity[:, axes] - velocity[:, axes]
            )))
            if (position_residual <= position_tolerance
                    and velocity_residual <= velocity_tolerance):
                return MostCoupledResponse(
                    platform_response, rotor_response, iteration,
                    position_residual, velocity_residual,
                )
            position = platform_response.position
            velocity = platform_response.velocity
        raise RuntimeError(
            "MOST platform/rotor coupling did not converge in "
            f"{max_iterations} passes (position {position_residual:.3g}, "
            f"velocity {velocity_residual:.3g})"
        )

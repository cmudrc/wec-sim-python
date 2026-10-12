"""Rotor and baseline-control dynamics for the published MOST turbine."""

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from wecsim.hydrodynamics.mostBEM import MostBEM
from wecsim.control.mostController import MostBaselineController
from wecsim.dynamics.mostMooring import MostStaticMooring
from wecsim.environment.turbSim import MostConstantWind, MostWindField


@dataclass(frozen=True)
class MostRotorResponse:
    time: np.ndarray
    rotor_speed: np.ndarray  # rad/s
    azimuth: np.ndarray  # rad
    generator_torque: np.ndarray  # N m
    blade_pitch: np.ndarray  # rad
    blade_root_load: np.ndarray  # N×6×3, preconed blade frames


@dataclass(frozen=True)
class MostRotor:
    """IEA 15 MW rotor on a prescribed six-DOF MOST platform trajectory.

    The platform pose and velocity are external inputs. Aerodynamic reaction
    and turbine mass do not feed back into platform motion in this runner.
    """

    bem: MostBEM
    controller: MostBaselineController
    inertia: float
    hub_offset: np.ndarray  # body-reference frame, m
    shaft_axis: np.ndarray  # body-reference unit vector
    initial_azimuth: float = 0.0

    @classmethod
    def from_iea15mw(cls, blade_directory: str | Path,
                     controller: MostBaselineController | None = None) -> "MostRotor":
        """Use the pinned MOST VolturnUS geometry and IEA 15 MW rotor."""
        if controller is None:
            controller = MostBaselineController.iea15mw()
        elif not isinstance(controller, MostBaselineController):
            raise TypeError("controller must be MostBaselineController")
        tilt = np.deg2rad(6)
        precone = np.deg2rad(4)
        blade_mass = 68507.602
        blade_cog = 27.601
        hub_radius = 3.94
        inertia = (1836784 + 973520 + 3 * np.cos(precone)**2 * (
            101086728 - blade_mass * blade_cog**2
            + blade_mass * (blade_cog + hub_radius)**2
        ))
        return cls(
            bem=MostBEM.from_iea15mw(blade_directory),
            controller=controller,
            inertia=float(inertia),
            hub_offset=np.array([
                -12.098 * np.cos(tilt), 0.0,
                15 + 129.386 + 4.3495 + 12.098 * np.sin(tilt) + 14.4,
            ]),
            shaft_axis=np.array([np.cos(tilt), 0.0, -np.sin(tilt)]),
        )

    def simulate(self, time, platform_position, platform_velocity,
                 wind: MostWindField | MostConstantWind) -> MostRotorResponse:
        """Advance rotor speed, azimuth, torque, and pitch from platform motion.

        ``platform_position`` and ``platform_velocity`` are N×6 world-frame
        body records in WEC-Sim's surge/sway/heave/roll/pitch/yaw order. The
        turbulent-wind input holds its final frame after the record ends.
        """
        time = np.asarray(time, dtype=float)
        position = np.asarray(platform_position, dtype=float)
        velocity = np.asarray(platform_velocity, dtype=float)
        if (time.ndim != 1 or time.size < 2
                or position.shape != (time.size, 6)
                or velocity.shape != (time.size, 6)
                or not np.isfinite(time).all()
                or not np.isfinite(position).all()
                or not np.isfinite(velocity).all()
                or not np.all(np.diff(time) > 0)):
            raise ValueError("MOST rotor needs aligned finite time and six-DOF platform records")
        if (not isinstance(wind, (MostWindField, MostConstantWind))
                or not np.isfinite(self.inertia) or self.inertia <= 0
                or np.shape(self.hub_offset) != (3,)
                or np.shape(self.shaft_axis) != (3,)):
            raise ValueError("invalid MOST rotor geometry or wind field")

        rotation = np.array([
            MostStaticMooring._rotation(*pose[3:]) for pose in position
        ])
        arm = rotation @ self.hub_offset
        shaft = rotation @ self.shaft_axis
        platform_shaft_speed = np.einsum("ij,ij->i", velocity[:, 3:], shaft)
        platform_shaft_acceleration = np.gradient(platform_shaft_speed, time)
        speed = np.empty(time.size)
        azimuth = np.empty(time.size)
        torque = np.empty(time.size)
        pitch = np.empty(time.size)
        blade_load = np.empty((time.size, 6, 3))
        speed[0] = self.controller.initial_omega
        azimuth[0] = self.initial_azimuth
        state = np.array([0.0, speed[0]/5, 0.0])
        torque[0], pitch[0] = self.controller._command(
            state, self.controller.initial_pitch,
        )
        precone = np.deg2rad(4)
        cosine, sine = np.cos(precone), np.sin(precone)

        def aerodynamic_torque(index, rotor_speed, rotor_azimuth, blade_pitch):
            hub_state = np.r_[
                position[index, :3] + arm[index], position[index, 3:],
                velocity[index, :3] + np.cross(velocity[index, 3:], arm[index]),
                velocity[index, 3:], rotor_azimuth, rotor_speed,
            ]
            loads = self.bem.loads(
                hub_state, blade_pitch, wind.sampler_at(float(time[index])),
            )
            if not np.isfinite(loads).all():
                raise ValueError("MOST BEM loads are not finite at the platform state")
            # Each BEM root load is in a preconed blade frame. Transfer its
            # moment to the shaft through the hub radius before projection.
            shaft_torque = np.sum(cosine*loads[3] - sine*loads[5]
                                  - self.bem.hub_radius*cosine*loads[1])
            return shaft_torque, loads

        for index in range(time.size - 1):
            step = time[index + 1] - time[index]
            aero, blade_load[index] = aerodynamic_torque(
                index, speed[index], azimuth[index], pitch[index],
            )
            acceleration = ((aero - torque[index]) / self.inertia
                            - platform_shaft_acceleration[index])
            predicted_speed = speed[index] + step*acceleration
            predicted_azimuth = azimuth[index] + step*speed[index]
            aero_next, _ = aerodynamic_torque(
                index + 1, predicted_speed, predicted_azimuth, pitch[index],
            )
            acceleration_next = ((aero_next - torque[index]) / self.inertia
                                 - platform_shaft_acceleration[index + 1])
            speed[index + 1] = speed[index] + step*(acceleration + acceleration_next)/2
            azimuth[index + 1] = azimuth[index] + step*(speed[index] + speed[index + 1])/2
            state, torque[index + 1], pitch[index + 1] = self.controller._advance(
                state, speed[index], speed[index + 1], step,
                torque[index], pitch[index],
            )
        _, blade_load[-1] = aerodynamic_torque(
            time.size - 1, speed[-1], azimuth[-1], pitch[-1],
        )
        return MostRotorResponse(time, speed, azimuth, torque, pitch,
                                 blade_load)

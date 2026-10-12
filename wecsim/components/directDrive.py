"""Coupled heave and simple direct-drive PTO dynamics.

The generator torque follows the requested reactive-controller force through
the winding L/R time constant. Gearbox inertia and friction act on the shaft,
and their reflected forces enter the body's heave equation.
"""

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp


@dataclass(frozen=True)
class DirectDriveResponse:
    time: np.ndarray
    position: np.ndarray
    velocity: np.ndarray
    acceleration: np.ndarray
    controller_force: np.ndarray
    shaft_velocity: np.ndarray
    shaft_torque: np.ndarray
    inertia_torque: np.ndarray
    friction_torque: np.ndarray
    generator_torque: np.ndarray
    current: np.ndarray
    voltage: np.ndarray
    resistance_loss: np.ndarray
    electrical_power: np.ndarray
    mechanical_power: np.ndarray
    body_force: np.ndarray


def integrate_direct_drive_heave(
    *, dt: float, end_time: float, height: float, period: float,
    ramp_time: float, rho: float, g: float, mass: float,
    displaced_volume: float, center_z: float, added_mass: float,
    radiation_damping: float, hydrostatic_stiffness: float,
    excitation_real: float, excitation_imaginary: float,
    kp: float, ki: float, torque_constant: float, gear_ratio: float,
    drivetrain_inertia: float, drivetrain_friction: float,
    winding_resistance: float, winding_inductance: float,
    initial_position: float | None = None, initial_velocity: float = 0.0,
) -> DirectDriveResponse:
    """Integrate one regular-wave heave body with a reactive direct drive."""
    values = [dt, end_time, height, period, ramp_time, rho, g, mass,
              displaced_volume, center_z, added_mass, radiation_damping,
              hydrostatic_stiffness, excitation_real, excitation_imaginary,
              kp, ki, torque_constant, gear_ratio, drivetrain_inertia,
              drivetrain_friction, winding_resistance, winding_inductance,
              initial_velocity]
    steps = round(end_time / dt)
    if (not np.isfinite(values).all() or dt <= 0 or end_time <= 0
            or height < 0 or period <= 0 or ramp_time < 0 or rho <= 0
            or g <= 0 or mass <= 0 or displaced_volume <= 0
            or torque_constant <= 0 or gear_ratio <= 0
            or drivetrain_inertia < 0 or drivetrain_friction < 0
            or winding_resistance <= 0 or winding_inductance <= 0
            or mass + added_mass + drivetrain_inertia * gear_ratio**2 <= 0
            or not np.isclose(steps * dt, end_time, rtol=0, atol=1e-10)):
        raise ValueError("direct-drive heave settings must be finite and physically valid")
    if initial_position is None:
        initial_position = center_z
    if not np.isfinite(initial_position):
        raise ValueError("initial position must be finite")

    omega = 2 * np.pi / period
    time_constant = winding_inductance / winding_resistance
    effective_mass = mass + added_mass + drivetrain_inertia * gear_ratio**2

    def derivative(t, state):
        z, speed, generator_torque = state
        ramp = (1.0 if ramp_time == 0 or t >= ramp_time else
                (1 - np.cos(np.pi * t / ramp_time)) / 2)
        excitation = height / 2 * ramp * (
            excitation_real * np.cos(omega * t)
            - excitation_imaginary * np.sin(omega * t)
        )
        controller = -kp * speed - ki * (z - center_z)
        acceleration = (
            excitation
            - hydrostatic_stiffness * (z - center_z)
            - (mass - rho * displaced_volume) * g
            - (radiation_damping + drivetrain_friction * gear_ratio**2) * speed
            + gear_ratio * generator_torque
        ) / effective_mass
        torque_rate = (controller / gear_ratio - generator_torque) / time_constant
        return [speed, acceleration, torque_rate]

    time = np.arange(steps + 1) * dt
    solution = solve_ivp(
        derivative, (0, end_time),
        [initial_position, initial_velocity, 0.0],
        t_eval=time, rtol=1e-8, atol=1e-10, max_step=dt,
    )
    if not solution.success or solution.y.shape != (3, len(time)):
        raise ValueError(f"direct-drive integration failed: {solution.message}")
    position, velocity, generator_torque = solution.y
    acceleration = np.array([
        derivative(t, solution.y[:, i])[1]
        for i, t in enumerate(time)
    ])
    controller_force = -kp * velocity - ki * (position - center_z)
    shaft_velocity = gear_ratio * velocity
    inertia_torque = -drivetrain_inertia * gear_ratio * acceleration
    friction_torque = -drivetrain_friction * shaft_velocity
    shaft_torque = inertia_torque + friction_torque + generator_torque
    current = generator_torque / torque_constant
    voltage = torque_constant * shaft_velocity
    resistance_loss = winding_resistance * current**2
    # The pinned Simulink block reports VI + I²R in this sign convention.
    electrical_power = voltage * current + resistance_loss
    mechanical_power = shaft_velocity * shaft_torque
    return DirectDriveResponse(
        time, position, velocity, acceleration, controller_force,
        shaft_velocity, shaft_torque, inertia_torque, friction_torque,
        generator_torque, current, voltage, resistance_loss,
        electrical_power, mechanical_power, gear_ratio * shaft_torque,
    )

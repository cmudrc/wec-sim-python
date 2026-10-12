"""Piecewise hydrodynamics for a single, vertically moving body.

The mass, equilibrium draft, and linear BEM coefficients change together at
specified times. This is the active motion in the published Variable_Mass
sphere application. The model keeps velocity continuous across a switch and
uses the combined rigid plus added mass in each state.
"""

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class HeaveHydroState:
    mass: float
    displaced_volume: float
    equilibrium_z: float
    added_mass: float
    radiation_damping: float
    hydrostatic_stiffness: float
    excitation_real: float
    excitation_imaginary: float


@dataclass(frozen=True)
class VariableHeaveResponse:
    time: np.ndarray
    position: np.ndarray
    velocity: np.ndarray
    acceleration: np.ndarray
    active_state: np.ndarray
    excitation: np.ndarray
    radiation: np.ndarray
    restoring: np.ndarray
    pto_force: np.ndarray


def integrate_variable_heave(
    states: tuple[HeaveHydroState, ...], switch_times: tuple[float, ...], *,
    dt: float, end_time: float, height: float, period: float, ramp_time: float,
    rho: float, g: float, pto_stiffness: float, pto_damping: float,
    pto_bias: float = 0.0, initial_position: float | None = None,
    initial_velocity: float = 0.0,
) -> VariableHeaveResponse:
    """Integrate a regular-wave heave model with complete state changes."""
    if len(states) < 2 or len(switch_times) != len(states) - 1:
        raise ValueError("variable heave needs one switch time between each state")
    switches = np.asarray(switch_times, dtype=float)
    steps = round(end_time / dt)
    if (not np.isfinite([dt, end_time, height, period, ramp_time, rho, g,
                         pto_stiffness, pto_damping, pto_bias,
                         initial_velocity]).all()
            or dt <= 0 or end_time <= 0 or period <= 0 or height < 0
            or ramp_time < 0 or rho <= 0 or g <= 0
            or not np.isclose(steps * dt, end_time, atol=1e-10)
            or switches.shape != (len(states) - 1,)
            or not np.isfinite(switches).all()
            or np.any(switches <= 0) or np.any(switches >= end_time)
            or np.any(np.diff(switches) <= 0)
            or not np.allclose(switches / dt, np.round(switches / dt),
                               rtol=0, atol=1e-8)):
        raise ValueError("variable heave settings must be finite and switches on the time grid")
    for state in states:
        values = tuple(vars(state).values())
        if (not np.isfinite(values).all() or state.mass <= 0
                or state.displaced_volume <= 0
                or state.mass + state.added_mass <= 0
                or state.radiation_damping < 0
                or state.hydrostatic_stiffness < 0):
            raise ValueError("each hydro state needs finite, physical heave coefficients")
    if initial_position is None:
        initial_position = states[0].equilibrium_z
    if not np.isfinite(initial_position):
        raise ValueError("initial position must be finite")

    omega = 2 * np.pi / period
    time = np.arange(steps + 1) * dt
    position = np.empty(steps + 1)
    velocity = np.empty(steps + 1)
    acceleration = np.empty(steps + 1)
    active_state = np.empty(steps + 1, dtype=int)
    excitation = np.empty(steps + 1)
    radiation = np.empty(steps + 1)
    restoring = np.empty(steps + 1)
    pto_force = np.empty(steps + 1)
    position[0] = initial_position
    velocity[0] = initial_velocity

    def forces(at_time, z, speed):
        # The tolerance resolves a floating-point value immediately below an
        # exact grid switch without moving the switch to an earlier step.
        index = int(np.searchsorted(switches, at_time + dt * 1e-8,
                                    side="right"))
        state = states[index]
        ramp = (1.0 if ramp_time == 0 or at_time >= ramp_time else
                (1 - np.cos(np.pi * at_time / ramp_time)) / 2)
        wave_force = height / 2 * ramp * (
            state.excitation_real * np.cos(omega * at_time)
            - state.excitation_imaginary * np.sin(omega * at_time)
        )
        damping_force = state.radiation_damping * speed
        restoring_force = (state.hydrostatic_stiffness
                           * (z - state.equilibrium_z)
                           + (state.mass - rho * state.displaced_volume) * g)
        pto = -pto_stiffness * (z - states[0].equilibrium_z) - pto_damping * speed + pto_bias
        a = ((wave_force - damping_force - restoring_force + pto)
             / (state.mass + state.added_mass))
        return index, a, wave_force, damping_force, restoring_force, pto

    for step, t in enumerate(time):
        (active_state[step], acceleration[step], excitation[step],
         radiation[step], restoring[step], pto_force[step]) = forces(
             t, position[step], velocity[step]
         )
        if step == steps:
            break

        def derivative(stage_time, z, speed):
            return np.array([speed, forces(stage_time, z, speed)[1]])

        z, speed = position[step], velocity[step]
        k1 = derivative(t, z, speed)
        k2 = derivative(t + dt/2, z + dt*k1[0]/2, speed + dt*k1[1]/2)
        k3 = derivative(t + dt/2, z + dt*k2[0]/2, speed + dt*k2[1]/2)
        k4 = derivative(t + dt, z + dt*k3[0], speed + dt*k3[1])
        position[step+1] = z + dt*(k1[0] + 2*k2[0] + 2*k3[0] + k4[0])/6
        velocity[step+1] = speed + dt*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1])/6

    return VariableHeaveResponse(time, position, velocity, acceleration,
                                 active_state, excitation, radiation,
                                 restoring, pto_force)

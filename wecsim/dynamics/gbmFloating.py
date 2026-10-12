"""Regular-wave floating body with rigid surge/heave/pitch and flexible modes.

The generalized mass is solved implicitly. MATLAB's sampled acceleration
feedback is a source-numerics diagnostic, not the default physical model.
"""

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from wecsim.hydrodynamics.bodyClass import BodyClass
from wecsim.environment.irregularWave import (
    IrregularComponents, synthesize_irregular_response,
)
from wecsim.components.orifice import OrificePTO


@dataclass(frozen=True)
class FloatingGBMResponse:
    time: np.ndarray
    body_position: np.ndarray
    body_velocity: np.ndarray
    mode_position: np.ndarray
    mode_velocity: np.ndarray
    mode_acceleration: np.ndarray
    wave_elevation: np.ndarray
    orifice_force: np.ndarray | None = None
    orifice_power: np.ndarray | None = None
    orifice_compressibility_flag: np.ndarray | None = None
    unwrapped_pitch: np.ndarray | None = None


def solve_floating_gbm_regular(
    hydro_file: str | Path, *, dt: float, end_time: float,
    height: float, period: float, pitch_inertia: float,
    mass: str | float = "equilibrium", ramp_time: float = 100,
    rho: float = 1000, g: float = 9.81,
) -> FloatingGBMResponse:
    """Solve one body with a 3-DOF floating joint and HDF5 flexible modes.

    The joint is at the body's origin; the equilibrium center of gravity must
    coincide with it. No PTO, mooring, drag, or body-to-body coupling is used.
    """
    body = BodyClass(str(hydro_file))
    body.bodyNumber = body.bodyTotal = 1
    body.readH5file()
    dof = int(np.asarray(body.dof).item())
    if dof <= 6 or int(np.asarray(body.dof_gbm).item()) != dof - 6:
        raise ValueError("floating GBM needs an HDF5 body with flexible modes")
    if not np.allclose(np.asarray(body.cg).ravel(), 0, atol=1e-10):
        raise ValueError("floating GBM currently needs the body CG at the joint origin")
    if mass == "equilibrium":
        physical_mass = rho * float(np.asarray(body.dispVol).item())
    else:
        physical_mass = float(mass)
    if physical_mass <= 0 or pitch_inertia <= 0:
        raise ValueError("body mass and pitch inertia must be positive")
    body.mass = mass
    omega = 2 * np.pi / period
    body.hydroForcePre(
        omega, [0], 1, np.array([0.0]), [], dt, rho, g,
        "regular", np.zeros((2, 1)), 1, 1, 0, 0, 0,
    )
    force = body.hydroForce
    active = np.r_[0, 2, 4, np.arange(6, dof)]
    n = len(active)
    mass_matrix = np.asarray(force["fAddedMass"])[np.ix_(active, active)].copy()
    mass_matrix[:3, :3] += np.diag([physical_mass, physical_mass, pitch_inertia])
    mass_matrix[3:, 3:] += force["gbm"]["mass_ff"]
    damping = np.asarray(force["fDamping"])[np.ix_(active, active)].copy()
    damping[3:, 3:] += force["gbm"]["damping"]
    stiffness = np.asarray(force["linearHydroRestCoef"])[np.ix_(active, active)].copy()
    stiffness[3:, 3:] += force["gbm"]["stiffness"]
    real = np.asarray(force["fExt"]["re"])[active]
    imaginary = np.asarray(force["fExt"]["im"])[active]
    inverse_mass = np.linalg.solve(mass_matrix, np.eye(n))

    def excitation(at_time):
        ramp = (1 if ramp_time == 0 or at_time >= ramp_time else
                (1 - np.cos(np.pi * at_time / ramp_time)) / 2)
        return height / 2 * ramp * (
            real * np.cos(omega * at_time)
            - imaginary * np.sin(omega * at_time)
        )

    def derivative(at_time, state):
        q, speed = state[:n], state[n:]
        acceleration = inverse_mass @ (
            excitation(at_time) - damping @ speed - stiffness @ q
        )
        return np.r_[speed, acceleration]

    time = np.arange(round(end_time / dt) + 1) * dt
    if not np.isclose(time[-1], end_time, rtol=0, atol=1e-10):
        raise ValueError("end_time must be an integer multiple of dt")
    state = np.zeros((len(time), 2 * n))
    for index in range(len(time) - 1):
        start = time[index]
        value = state[index]
        k1 = derivative(start, value)
        k2 = derivative(start + dt / 2, value + dt * k1 / 2)
        k3 = derivative(start + dt / 2, value + dt * k2 / 2)
        k4 = derivative(start + dt, value + dt * k3)
        state[index + 1] = value + dt * (k1 + 2*k2 + 2*k3 + k4) / 6
    acceleration = np.array([
        derivative(at_time, value)[n:]
        for at_time, value in zip(time, state)
    ])
    body_position = np.zeros((len(time), 6))
    body_velocity = np.zeros_like(body_position)
    body_position[:, [0, 2, 4]] = state[:, :3]
    body_velocity[:, [0, 2, 4]] = state[:, n:n+3]
    ramp = np.ones_like(time)
    if ramp_time:
        early = time < ramp_time
        ramp[early] = (1 - np.cos(np.pi * time[early] / ramp_time)) / 2
    return FloatingGBMResponse(
        time, body_position, body_velocity,
        state[:, 3:n], state[:, n+3:], acceleration[:, 3:],
        height / 2 * ramp * np.cos(omega * time),
    )


def solve_floating_gbm_pm_orifice(
    hydro_file: str | Path, *, dt: float, end_time: float,
    components: IrregularComponents, pitch_inertia: float,
    orifice: OrificePTO, mass: str | float = "equilibrium",
    orifice_force_path: str = "coupled",
    ramp_time: float = 10, memory_time: float = 15,
    heave_linear_damping: float = 0, mode_linear_damping: float = 0,
    heave_drag_cd: float = 0, heave_drag_area: float = 0,
    pitch_drag_cd: float = 0, pitch_drag_area: float = 0,
    rho: float = 1000, g: float = 9.81,
) -> FloatingGBMResponse:
    """Couple one OWC flexible mode and its orifice to a floating body.

    The rigid joint permits surge, heave, and pitch. The default applies equal
    and opposite piston forces to the flexible state and rigid heave.
    ``published_owc`` reproduces the pinned Simulink signal route, where the
    piston affects reported flexible acceleration and rigid added-mass
    feedback but is absent from the flexible state-space input. It is for
    source comparisons, not physical predictions. Rigid pitch restoring uses
    the source's Euler pitch branch after 90 degrees. Returned flexible
    acceleration is always the derivative of the simulated state; the pinned
    source reports a separate piston-inclusive acceleration signal.
    """
    if not isinstance(orifice, OrificePTO):
        raise TypeError("orifice must be an OrificePTO")
    if (not isinstance(orifice_force_path, str)
            or orifice_force_path not in ("coupled", "published_owc")):
        raise ValueError("orifice_force_path must be coupled or published_owc")
    if (not np.isfinite([dt, end_time, ramp_time, memory_time,
                         pitch_inertia, rho, g,
                         heave_linear_damping, mode_linear_damping,
                         heave_drag_cd, heave_drag_area,
                         pitch_drag_cd, pitch_drag_area]).all()
            or dt <= 0 or end_time <= 0 or ramp_time < 0
            or memory_time <= 0 or pitch_inertia <= 0 or rho <= 0 or g <= 0
            or min(heave_linear_damping, mode_linear_damping,
                   heave_drag_cd, heave_drag_area,
                   pitch_drag_cd, pitch_drag_area) < 0):
        raise ValueError("floating OWC settings must be finite and nonnegative")
    if (not np.isclose(end_time / dt, round(end_time / dt), rtol=0, atol=1e-10)
            or not np.isclose(memory_time / dt, round(memory_time / dt),
                              rtol=0, atol=1e-10)):
        raise ValueError("end_time and memory_time must be multiples of dt")
    body = BodyClass(str(hydro_file))
    body.bodyNumber = body.bodyTotal = 1
    body.readH5file()
    if int(np.asarray(body.dof).item()) != 7 or int(np.asarray(body.dof_gbm).item()) != 1:
        raise ValueError("floating OWC needs one rigid body and one flexible mode")
    if not np.allclose(np.asarray(body.cg).ravel(), 0, atol=1e-10):
        raise ValueError("floating OWC needs the body CG at its joint origin")
    physical_mass = (rho * float(np.asarray(body.dispVol).item())
                     if mass == "equilibrium" else float(mass))
    if not np.isfinite(physical_mass) or physical_mass <= 0:
        raise ValueError("body mass must be positive and finite")
    body.mass = mass
    body.linearDamping = np.diag([
        0, 0, heave_linear_damping, 0, 0, 0, mode_linear_damping,
    ])
    body.viscDrag["cd"] = np.array([
        0, 0, heave_drag_cd, 0, pitch_drag_cd, 0, 0,
    ])
    body.viscDrag["characteristicArea"] = np.array([
        0, 0, heave_drag_area, 0, pitch_drag_area, 0, 0,
    ])
    incident = synthesize_irregular_response(
        hydro_file, components, dt=dt, end_time=end_time,
        ramp_time=ramp_time, rho=rho, g=g,
    )
    memory = np.arange(round(memory_time / dt) + 1) * dt
    body.hydroForcePre(
        components.omega, [0], len(memory), memory,
        len(components.omega), dt, rho, g, "irregular",
        np.zeros((2, 1)), 1, 1, 0, 0, 0,
    )
    force = body.hydroForce
    active = np.array([0, 2, 4, 6])
    added_mass = np.asarray(force["fAddedMass"])[np.ix_(active, active)]
    matrix = added_mass.copy()
    matrix[:3, :3] += np.diag([physical_mass, physical_mass, pitch_inertia])
    matrix[3, 3] = float(np.asarray(force["gbm"]["mass_ff"])[0, 0])
    if matrix[3, 3] <= 0:
        raise ValueError("floating OWC flexible effective mass must be positive")
    # In the source model, heave added mass sees the reported flexible
    # acceleration (state acceleration + piston force / flexible mass), even
    # though the flexible state itself receives no piston force.
    piston_distribution = (
        np.array([0, -(1 + matrix[1, 3] / matrix[3, 3]), 0, 0])
        if orifice_force_path == "published_owc"
        else np.array([0, -1, 0, 1])
    )
    stiffness = np.asarray(force["linearHydroRestCoef"])[np.ix_(active, active)].copy()
    stiffness[3, 3] += float(np.asarray(force["gbm"]["stiffness"])[0, 0])
    damping = np.asarray(force["linearDamping"])[np.ix_(active, active)].copy()
    damping[3, 3] += float(np.asarray(force["gbm"]["damping"])[0, 0])
    drag = np.diag(np.asarray(force["visDrag"])[np.ix_(active, active)])
    kernel = np.asarray(force["irkb"])[:, active, :][:, :, active]
    excitation = incident.excitation_force[:, active]
    count = len(incident.time)
    coordinate = np.zeros((count, 4))
    velocity = np.zeros_like(coordinate)
    acceleration = np.zeros_like(coordinate)
    pto_force = np.zeros(count)
    pto_power = np.zeros(count)
    pto_flag = np.zeros(count, dtype=bool)
    pto_coefficient = (0.5 * orifice.air_density * orifice.piston_area
                       * (orifice.piston_area /
                          (orifice.discharge_coefficient * orifice.orifice_area))**2)
    for step in range(1, count):
        lag = min(step, len(kernel) - 1)
        radiation = dt * np.einsum(
            "kij,kj->i", kernel[1:lag + 1],
            velocity[step - lag:step][::-1],
        )
        if step >= len(kernel) - 1:
            radiation -= dt / 2 * kernel[-1] @ velocity[step - lag]
        previous_q = coordinate[step - 1]
        previous_v = velocity[step - 1]
        previous_a = acceleration[step - 1]
        trial = previous_v.copy()
        for _ in range(12):
            position = previous_q + dt / 2 * (previous_v + trial)
            reaction = float(orifice.evaluate(trial[3]).force)
            applied = piston_distribution * reaction
            restoring = stiffness @ position
            restoring[2] = stiffness[2, 2] * np.arcsin(np.sin(position[2]))
            current = (excitation[step] - radiation
                       - dt / 2 * kernel[0] @ trial
                       - restoring - damping @ trial
                       - drag * trial * np.abs(trial) + applied)
            residual = (2 / dt * matrix @ (trial - previous_v)
                        - matrix @ previous_a - current)
            jacobian = (2 / dt * matrix + dt / 2 * (kernel[0] + stiffness)
                        + damping + np.diag(2 * drag * np.abs(trial)))
            jacobian[2, 2] += (dt / 2 * stiffness[2, 2]
                               * (np.sign(np.cos(position[2])) - 1))
            jacobian[:, 3] += piston_distribution * (
                2 * pto_coefficient * abs(trial[3])
            )
            correction = np.linalg.solve(jacobian, residual)
            trial -= correction
            if np.max(np.abs(correction)) < 1e-11:
                break
        else:
            raise RuntimeError(f"floating OWC step {step} did not converge")
        velocity[step] = trial
        coordinate[step] = previous_q + dt / 2 * (previous_v + trial)
        acceleration[step] = 2 / dt * (trial - previous_v) - previous_a
        result = orifice.evaluate(trial[3])
        pto_force[step] = result.force
        pto_power[step] = result.absorbed_power
        pto_flag[step] = result.compressibility_flag
    body_position = np.zeros((count, 6))
    body_velocity = np.zeros_like(body_position)
    body_position[:, [0, 2, 4]] = coordinate[:, :3]
    body_position[:, 4] = np.arcsin(np.sin(coordinate[:, 2]))
    inverted_branch = np.cos(coordinate[:, 2]) < 0
    body_position[inverted_branch, 3] = np.pi
    body_position[inverted_branch, 5] = np.pi
    body_velocity[:, [0, 2, 4]] = velocity[:, :3]
    return FloatingGBMResponse(
        incident.time, body_position, body_velocity,
        coordinate[:, 3:], velocity[:, 3:], acceleration[:, 3:],
        incident.elevation, pto_force, pto_power, pto_flag,
        coordinate[:, 2],
    )

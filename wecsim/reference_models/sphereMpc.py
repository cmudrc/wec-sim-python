"""Published single-body Sphere heave case with model predictive PTO control."""

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from wecsim.hydrodynamics.bodyClass import BodyClass
from wecsim.dynamics.generalDynamics import BodyMotion, DynamicBody, GeneralizedDynamics
from wecsim.environment.irregularWave import (
    jonswap_equal_energy_components, synthesize_irregular_response,
)
from wecsim.control.mpcControl import simulate_sphere_mpc
from wecsim.control.mpcPlant import build_sphere_mpc_matrices


@dataclass(frozen=True)
class SphereMPCResult:
    time: np.ndarray
    position: np.ndarray
    velocity: np.ndarray
    wave_elevation: np.ndarray
    excitation_force: np.ndarray
    pto_force: np.ndarray
    absorbed_power: np.ndarray
    command_rate: np.ndarray
    internal_state: np.ndarray
    feasible: np.ndarray


def run_sphere_mpc(
    hydro_file: str | Path, coefficient_file: str | Path, *,
    significant_height: float = 2.5, peak_period: float = 8,
    phase: np.ndarray | None = None, seed: int | None = None,
    gamma: float | None = None, frequency_count: int = 500,
    dt: float = .01, end_time: float = 400,
    ramp_time: float = 100, radiation_memory: float = 10,
    control_step: float = .5, prediction_horizon: float = 15,
    command_delay: float = .5, start_time: float = 205,
    history_steps: int = 200, order: int = 4,
    rate_penalty: float = 1e-7,
    max_speed: float = 3, max_position: float = 4,
    max_force: float = 2e6, max_force_rate: float = 1.5e6,
    rho: float = 1000, g: float = 9.81,
) -> SphereMPCResult:
    """Simulate the published heave-only Sphere and its MPC PTO.

    The source's prediction plant is driven by wave excitation independently
    of the physical body. Its force-rate command passes through a 0.5 s
    hold/transition before the PTO force integrator. Settings can be changed,
    but paired MATLAB trajectory gates cover the published defaults only.
    """
    components = jonswap_equal_energy_components(
        hydro_file, significant_height=significant_height,
        peak_period=peak_period, directions=np.array([0.0]),
        spreading=np.array([1.0]), count=frequency_count,
        phase=phase, seed=seed, gamma=gamma,
    )
    incident = synthesize_irregular_response(
        hydro_file, components, dt=dt, end_time=end_time,
        ramp_time=ramp_time, rho=rho, g=g,
    )
    matrices = build_sphere_mpc_matrices(
        hydro_file, coefficient_file, prediction_step=control_step,
        prediction_horizon=prediction_horizon, rate_penalty=rate_penalty,
        rho=rho, g=g,
    )
    controller = simulate_sphere_mpc(
        matrices, incident.excitation_force[:, 2], dt=dt,
        control_step=control_step, command_delay=command_delay,
        start_time=start_time, history_steps=history_steps, order=order,
        max_speed=max_speed, max_position=max_position,
        max_force=max_force, max_force_rate=max_force_rate,
    )

    body = BodyClass(str(hydro_file))
    body.bodyNumber = body.bodyTotal = 1
    body.readH5file()
    body.mass = "equilibrium"
    body.hydroStiffness = np.zeros((6, 6))
    body.viscDrag = {
        "Drag": np.zeros((6, 6)), "cd": np.zeros(6),
        "characteristicArea": np.zeros(6),
    }
    body.linearDamping = np.zeros((6, 6))
    memory_steps = round(radiation_memory / dt)
    if (dt <= 0 or radiation_memory <= 0
            or not np.isclose(memory_steps * dt, radiation_memory, atol=1e-10)):
        raise ValueError("radiation_memory must be a positive multiple of dt")
    memory_time = np.arange(memory_steps + 1) * dt
    wave_amplitude = np.vstack((incident.time, np.zeros_like(incident.time)))
    body.hydroForcePre(
        components.omega, components.directions,
        len(memory_time), memory_time, len(components.omega),
        dt, rho, g, "irregular", wave_amplitude, 1, 1, 0, 0, 0,
    )
    mapping = np.zeros((6, 1))
    mapping[2, 0] = 1

    def motion(coordinate, speed):
        displacement = np.zeros(6)
        displacement[2] = coordinate[0]
        return BodyMotion(displacement, mapping, np.zeros(6))

    rigid_mass = np.zeros((6, 6))
    rigid_mass[:3, :3] = np.eye(3) * float(np.asarray(body.mass).item())
    equilibrium = np.r_[np.asarray(body.cg).ravel(), 0, 0, 0]
    dynamic = DynamicBody(
        rigid_mass=rigid_mass,
        added_mass=(np.asarray(body.hydroForce["fAddedMass"]),),
        damping=(np.zeros((6, 6)),),
        restoring=np.asarray(body.hydroForce["linearHydroRestCoef"]),
        static_force=np.zeros(6), reference_position=equilibrium,
        motion=motion,
        excitation=lambda at: incident.excitation_force[round(at / dt)],
        radiation_kernel=np.asarray(body.hydroForce["irkb"]),
    )
    physical = GeneralizedDynamics((dynamic,), 1).integrate(
        dt=dt, end_time=end_time,
        applied_force_history=controller.pto_force[:, None],
    )
    velocity = physical.body_velocity[:, 0, 2]
    return SphereMPCResult(
        time=physical.time,
        position=physical.body_position[:, 0, 2],
        velocity=velocity,
        wave_elevation=incident.elevation,
        excitation_force=incident.excitation_force[:, 2],
        pto_force=controller.pto_force,
        absorbed_power=-controller.pto_force * velocity,
        command_rate=controller.command_rate,
        internal_state=controller.internal_state,
        feasible=controller.feasible,
    )

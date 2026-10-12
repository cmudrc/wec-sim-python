"""Coupled pitch and reverse-osmosis dynamics for published OSWEC Desalination."""

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

from wecsim.hydrodynamics.bodyClass import BodyClass
from wecsim.components.crank import PitchRodLinkage
from wecsim.components.desalination import FourValveRectifiedCylinder, ReverseOsmosisHydraulicNetwork
from wecsim.dynamics.generalDynamics import BodyMotion, DynamicBody, GeneralizedDynamics
from wecsim.environment.irregularWave import IrregularComponents, synthesize_irregular_response
from wecsim.hydrodynamics.morison import MorisonElement, irregular_morison_source_drag


@dataclass(frozen=True)
class OSWECDesalinationResponse:
    time: np.ndarray
    pitch: np.ndarray
    pitch_velocity: np.ndarray
    body_position: np.ndarray
    body_velocity: np.ndarray
    rod_speed: np.ndarray
    rod_force: np.ndarray
    chamber_pressure: np.ndarray
    high_pressure: np.ndarray
    permeate_flow: np.ndarray
    relief_flow: np.ndarray
    wave_elevation: np.ndarray
    excitation_force: np.ndarray


def run_oswec_desalination(
    hydro_file: str | Path,
    components: IrregularComponents,
    network: ReverseOsmosisHydraulicNetwork,
    valves: FourValveRectifiedCylinder,
    linkage: PitchRodLinkage,
    morison_elements: Sequence[MorisonElement],
    *,
    dt: float = .01,
    end_time: float = 300.,
    ramp_time: float = 50.,
    radiation_memory: float = 30.,
    water_depth: float = 10.9,
    body_mass: float = 127_000.,
    pitch_inertia: float = 1.85e6,
    hinge_z: float = -8.9,
    initial_pitch: float = 0.,
    initial_speed: float = 0.,
    rho: float = 1000.,
    g: float = 9.81,
) -> OSWECDesalinationResponse:
    """Advance flap, radiation memory, Morison drag, and hydraulic PTO.

    The published actuation block delays measured cylinder force by one
    sample. The four-valve chamber model uses passive orifice equations, so
    it does not reproduce the legacy source's alternating raw pressures.
    Saved MATLAB phases may be supplied as ``components.phase`` for a paired
    incident sea; all motion and hydraulic states are computed in Python.
    """
    if (not isinstance(components, IrregularComponents)
            or not isinstance(network, ReverseOsmosisHydraulicNetwork)
            or not isinstance(valves, FourValveRectifiedCylinder)
            or not isinstance(linkage, PitchRodLinkage)
            or not morison_elements):
        raise TypeError("desalination runner needs sea, network, valves, linkage, and Morison elements")
    if (not np.isclose(network.cylinder_area, valves.cylinder.area_a)
            or not np.isclose(network.cylinder_area, valves.cylinder.area_b)):
        raise ValueError("desalination rectifier needs equal cylinder and feed areas")
    settings = np.asarray((
        dt, end_time, ramp_time, radiation_memory, water_depth, body_mass,
        pitch_inertia, hinge_z, initial_pitch, initial_speed, rho, g,
    ), dtype=float)
    if (not np.isfinite(settings).all() or dt <= 0 or end_time <= 0
            or ramp_time < 0 or radiation_memory < dt or water_depth <= 0
            or body_mass <= 0 or pitch_inertia <= 0 or rho <= 0 or g <= 0):
        raise ValueError("desalination time, geometry, and body settings are invalid")
    steps = round(end_time / dt)
    memory_steps = round(radiation_memory / dt)
    if (not np.isclose(steps * dt, end_time, rtol=0, atol=1e-10)
            or not np.isclose(memory_steps * dt, radiation_memory,
                              rtol=0, atol=1e-10)):
        raise ValueError("duration and radiation memory must be integer steps")

    count = steps + 1
    time = np.arange(count) * dt
    sea = synthesize_irregular_response(
        hydro_file, components, dt=dt, end_time=end_time,
        ramp_time=ramp_time, body_number=1, rho=rho, g=g,
    )
    wave_force = sea.excitation_force[:, :6]

    body = BodyClass(str(hydro_file))
    body.bodyNumber = 1
    body.readH5file()
    if int(np.asarray(body.dof).item()) != 6:
        raise ValueError("desalination flap needs six hydrodynamic DOFs")
    cg_z = float(np.asarray(body.cg).ravel()[2])
    lever = cg_z - hinge_z
    if lever <= 0:
        raise ValueError("the hinge must be below the flap center of gravity")
    irf_time = body.hydroData["hydro_coeffs"]["radiation_damping"][
        "impulse_response_fun"]["t"]
    if radiation_memory > np.max(irf_time) + 1e-10:
        raise ValueError("radiation memory exceeds the HDF5 kernel")
    convolution_time = np.arange(memory_steps + 1) * dt
    body.mass = body_mass
    body.hydroStiffness = np.zeros((6, 6))
    body.viscDrag = {"Drag": np.zeros((6, 6)), "cd": np.zeros(6),
                     "characteristicArea": np.zeros(6)}
    body.linearDamping = np.zeros((6, 6))
    body.hydroForcePre(
        [], [0], len(convolution_time), convolution_time, [], dt, rho, g,
        "noWaveCIC", np.vstack((time, np.zeros_like(time))),
        1, 1, 0, 0, 0,
    )
    hydro = body.hydroForce
    kernel = np.asarray(hydro["irkb"])

    def motion(q, v):
        angle, speed = q[0], v[0]
        jacobian = np.array([
            lever * np.cos(angle), 0., -lever * np.sin(angle),
            0., 1., 0.,
        ])[:, None]
        bias = np.array([
            -lever * np.sin(angle), 0., -lever * np.cos(angle),
            0., 0., 0.,
        ]) * speed**2
        displacement = np.array([
            lever * np.sin(angle), 0., lever * (np.cos(angle) - 1.),
            0., angle, 0.,
        ])
        return BodyMotion(displacement, jacobian, bias)

    reference = np.array([0., 0., cg_z, 0., 0., 0.])

    def state_excitation(at_time, q, v):
        state = motion(q, v)
        return (wave_force[round(at_time / dt)]
                + irregular_morison_source_drag(
                    morison_elements, time=at_time,
                    position=reference + state.displacement,
                    velocity=state.jacobian[:, 0] * v[0],
                    components=components, water_depth=water_depth,
                    ramp_time=ramp_time, rho=rho, g=g,
                ))

    device = GeneralizedDynamics((DynamicBody(
        rigid_mass=np.diag([body_mass, body_mass, body_mass,
                            0., pitch_inertia, 0.]),
        added_mass=(np.asarray(hydro["fAddedMass"]),),
        damping=(np.zeros((6, 6)),),
        restoring=np.asarray(hydro["linearHydroRestCoef"]),
        static_force=np.array([
            0., 0., (rho * float(np.asarray(body.dispVol).item())
                       - body_mass) * g, 0., 0., 0.,
        ]),
        reference_position=reference,
        motion=motion, excitation=lambda at_time: wave_force[round(at_time / dt)],
        state_excitation=state_excitation, radiation_kernel=kernel,
    ),), 1)

    angle = np.zeros(count)
    speed = np.zeros(count)
    acceleration = np.zeros(count)
    velocity_history = np.zeros((count, 6))
    cylinder_speed = np.empty(count)
    rod_force = np.empty(count)
    chamber_pressure = np.empty((count, 2))
    high_pressure = np.empty(count)
    permeate_flow = np.empty(count)
    relief_flow = np.empty(count)
    angle[0], speed[0] = initial_pitch, initial_speed
    hydraulic = network.initial_state()

    def record_hydraulic(step):
        jacobian = linkage.stroke_and_jacobian(angle[step])[1]
        cylinder_speed[step] = -jacobian * speed[step]
        high_pressure[step] = hydraulic.pressure
        permeate_flow[step] = hydraulic.permeate_flow
        relief_flow[step] = hydraulic.relief_flow
        chamber_pressure[step] = valves.chamber_pressures(
            cylinder_speed[step], hydraulic.pressure,
        )
        rod_force[step] = valves.cylinder.force(*chamber_pressure[step])

    record_hydraulic(0)
    acceleration[0] = device.acceleration(
        time[0], np.array([angle[0]]), np.array([speed[0]]),
        known_radiation=(np.zeros(6),), dt=dt,
        applied_force=np.array([0.]),
    )[0]
    velocity_history[0] = (
        motion(np.array([angle[0]]), np.array([speed[0]])).jacobian[:, 0]
        * speed[0]
    )

    for step in range(1, count):
        memory = min(step, memory_steps)
        known = dt * np.einsum(
            "tij,tj->i", kernel[1:memory + 1],
            velocity_history[step - memory:step][::-1],
        )
        if step >= memory_steps:
            known -= dt / 2 * kernel[-1] @ velocity_history[step - memory]
        trial_speed = speed[step - 1]
        for _ in range(30):
            trial_angle = (angle[step - 1]
                           + dt * (speed[step - 1] + trial_speed) / 2)
            jacobian = linkage.stroke_and_jacobian(trial_angle)[1]
            torque = -rod_force[step - 1] * jacobian
            trial_acceleration = device.acceleration(
                time[step], np.array([trial_angle]),
                np.array([trial_speed]), known_radiation=(known,), dt=dt,
                applied_force=np.array([torque]),
            )[0]
            next_speed = (speed[step - 1]
                          + dt * (acceleration[step - 1]
                                  + trial_acceleration) / 2)
            if abs(next_speed - trial_speed) < 1e-12:
                trial_speed = next_speed
                break
            trial_speed = next_speed
        else:
            raise RuntimeError("desalination radiation-memory step did not converge")
        speed[step] = trial_speed
        angle[step] = (angle[step - 1]
                       + dt * (speed[step - 1] + speed[step]) / 2)
        jacobian = linkage.stroke_and_jacobian(angle[step])[1]
        acceleration[step] = device.acceleration(
            time[step], np.array([angle[step]]), np.array([speed[step]]),
            known_radiation=(known,), dt=dt,
            applied_force=np.array([-rod_force[step - 1] * jacobian]),
        )[0]
        velocity_history[step] = (
            motion(np.array([angle[step]]), np.array([speed[step]])).jacobian[:, 0]
            * speed[step]
        )
        hydraulic = network.step(-jacobian * speed[step], hydraulic, dt)
        record_hydraulic(step)

    center_position = np.column_stack((
        lever * np.sin(angle), np.zeros(count),
        cg_z + lever * (np.cos(angle) - 1),
    ))
    return OSWECDesalinationResponse(
        time=time, pitch=angle, pitch_velocity=speed,
        body_position=center_position, body_velocity=velocity_history[:, :3],
        rod_speed=cylinder_speed, rod_force=rod_force,
        chamber_pressure=chamber_pressure, high_pressure=high_pressure,
        permeate_flow=permeate_flow, relief_flow=relief_flow,
        wave_elevation=sea.elevation, excitation_force=wave_force,
    )

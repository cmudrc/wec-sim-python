"""Coupled OSWEC pitch and rectified hydraulic PTO for two published layouts."""

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from wecsim.hydrodynamics.bodyClass import BodyClass
from wecsim.components.crank import AdjustableRodCrank, FixedRodCrank
from wecsim.dynamics.generalDynamics import BodyMotion, DynamicBody, GeneralizedDynamics
from wecsim.components.hydraulic import RectifiedHydraulicPTO
from wecsim.environment.irregularWave import IrregularComponents, synthesize_irregular_response


@dataclass(frozen=True)
class OSWECHydraulicResponse:
    time: np.ndarray
    pitch: np.ndarray
    pitch_velocity: np.ndarray
    body_position: np.ndarray  # time, xyz at the flap center
    body_velocity: np.ndarray
    crank_stroke: np.ndarray
    crank_speed: np.ndarray
    crank_torque: np.ndarray
    cylinder_force: np.ndarray
    cylinder_pressure: np.ndarray  # time, chamber A/B
    accumulator_pressure: np.ndarray  # time, high/low
    valve_flow: np.ndarray  # time, ports A/B/high/low
    motor_torque: np.ndarray
    motor_flow: np.ndarray
    shaft_speed_rpm: np.ndarray
    generator_current: np.ndarray
    generator_voltage: np.ndarray
    load_resistance: np.ndarray
    wave_elevation: np.ndarray
    excitation_force: np.ndarray  # time, six flap axes


def run_oswec_rectified_hydraulic(
    hydro_file: str | Path,
    components: IrregularComponents,
    pto: RectifiedHydraulicPTO,
    linkage: FixedRodCrank | AdjustableRodCrank,
    *,
    dt: float = .01,
    end_time: float = 400.,
    ramp_time: float = 100.,
    radiation_memory: float = 30.,
    body_mass: float = 127_000.,
    pitch_inertia: float = 1.85e6,
    pitch_damping: float = 1e7,
    hinge_z: float = -8.9,
    initial_angle: float = 0.,
    initial_speed: float = 0.,
    rho: float = 1000.,
    g: float = 9.81,
) -> OSWECHydraulicResponse:
    """Integrate the published fixed-base OSWEC hydraulic PTO layouts.

    The linkage converts pitch to cylinder stroke and force to hinge torque.
    The seven hydraulic/electrical states use the source's forward Euler
    updates; flap pitch uses an implicit trapezoidal step with radiation IRF.
    Saved MATLAB phases can be supplied for an exact incident-wave pair.
    """
    if (not isinstance(components, IrregularComponents)
            or not isinstance(pto, RectifiedHydraulicPTO)
            or not isinstance(linkage, (FixedRodCrank, AdjustableRodCrank))):
        raise TypeError("OSWEC hydraulic runner needs sea, PTO, and crank objects")
    values = np.asarray((dt, end_time, ramp_time, radiation_memory, body_mass,
                         pitch_inertia, pitch_damping, hinge_z, initial_angle,
                         initial_speed, rho, g), dtype=float)
    if (not np.isfinite(values).all() or dt <= 0 or end_time <= 0
            or ramp_time < 0 or radiation_memory <= 0 or body_mass <= 0
            or pitch_inertia <= 0 or pitch_damping < 0 or rho <= 0 or g <= 0):
        raise ValueError("OSWEC hydraulic time and body settings are invalid")
    steps = round(end_time / dt)
    memory_steps = round(radiation_memory / dt)
    if (not np.isclose(steps * dt, end_time, rtol=0, atol=1e-10)
            or not np.isclose(memory_steps * dt, radiation_memory,
                              rtol=0, atol=1e-10)):
        raise ValueError("duration and radiation memory must be integer steps")
    time = np.arange(steps + 1) * dt
    sea = synthesize_irregular_response(
        hydro_file, components, dt=dt, end_time=end_time,
        ramp_time=ramp_time, body_number=1, rho=rho, g=g,
    )
    excitation_force = sea.excitation_force[:, :6]

    body = BodyClass(str(hydro_file))
    body.bodyNumber = 1
    body.readH5file()
    if int(np.asarray(body.dof).item()) != 6:
        raise ValueError("OSWEC hydraulic runner needs a six-DOF flap")
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
    body.linearDamping = np.diag([0, 0, 0, 0, pitch_damping, 0])
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

    def excitation(at_time):
        return excitation_force[round(at_time / dt)]

    device = GeneralizedDynamics((DynamicBody(
        rigid_mass=np.diag([body_mass, body_mass, body_mass,
                            0., pitch_inertia, 0.]),
        added_mass=(np.asarray(hydro["fAddedMass"]),),
        damping=(np.asarray(body.linearDamping),),
        restoring=np.asarray(hydro["linearHydroRestCoef"]),
        static_force=np.array([
            0., 0., (rho * float(np.asarray(body.dispVol).item())
                       - body_mass) * g, 0., 0., 0.,
        ]),
        reference_position=np.array([0., 0., cg_z, 0., 0., 0.]),
        motion=motion, excitation=excitation, radiation_kernel=kernel,
    ),), 1)

    angle = np.zeros(steps + 1)
    speed = np.zeros(steps + 1)
    acceleration = np.zeros(steps + 1)
    angle[0], speed[0] = initial_angle, initial_speed
    states = np.empty((steps + 1, 7))
    states[0] = pto.initial_state()
    cylinder_force = np.empty(steps + 1)
    cylinder_force[0] = pto.force(states[0])
    torque = np.empty(steps + 1)
    torque[0] = linkage.torque(cylinder_force[0], angle[0])
    velocity_history = np.zeros((steps + 1, 6))
    zeros = (np.zeros(6),)
    acceleration[0] = device.acceleration(
        time[0], angle[:1], speed[:1], known_radiation=zeros, dt=dt,
        applied_force=np.array([torque[0]]),
    )[0]
    velocity_history[0] = motion(angle[:1], speed[:1]).jacobian[:, 0] * speed[0]

    for step in range(1, steps + 1):
        previous = step - 1
        stroke, jacobian = linkage.stroke_and_jacobian(angle[previous])
        rates = pto.state_rate(
            stroke, jacobian * speed[previous], states[previous],
        )
        states[step] = states[previous] + dt * rates
        if not np.isfinite(states[step]).all():
            raise FloatingPointError("OSWEC hydraulic PTO state became nonfinite")
        cylinder_force[step] = pto.force(states[step])
        memory = min(step, memory_steps)
        known = dt * np.einsum(
            "tij,tj->i", kernel[1:memory + 1],
            velocity_history[step - memory:step][::-1],
        )
        if step >= memory_steps:
            known -= dt / 2 * kernel[-1] @ velocity_history[step - memory]
        trial_speed = speed[previous]
        for _ in range(30):
            trial_angle = angle[previous] + dt * (speed[previous] + trial_speed) / 2
            trial_torque = linkage.torque(cylinder_force[step], trial_angle)
            trial_acceleration = device.acceleration(
                time[step], np.array([trial_angle]), np.array([trial_speed]),
                known_radiation=(known,), dt=dt,
                applied_force=np.array([trial_torque]),
            )[0]
            next_speed = speed[previous] + dt * (
                acceleration[previous] + trial_acceleration
            ) / 2
            if abs(next_speed - trial_speed) < 1e-12:
                trial_speed = next_speed
                break
            trial_speed = next_speed
        else:
            raise RuntimeError("OSWEC hydraulic step did not converge; reduce dt")
        speed[step] = trial_speed
        angle[step] = angle[previous] + dt * (speed[previous] + speed[step]) / 2
        torque[step] = linkage.torque(cylinder_force[step], angle[step])
        acceleration[step] = device.acceleration(
            time[step], np.array([angle[step]]), np.array([speed[step]]),
            known_radiation=(known,), dt=dt,
            applied_force=np.array([torque[step]]),
        )[0]
        velocity_history[step] = (
            motion(np.array([angle[step]]), np.array([speed[step]])).jacobian[:, 0]
            * speed[step]
        )

    crank_stroke, crank_jacobian = linkage.stroke_and_jacobian(angle)
    pressure_high = pto.high_accumulator.pressure(states[:, 2])
    pressure_low = pto.low_accumulator.pressure(states[:, 3])
    valve_flow = np.column_stack(pto.valve.flows(
        states[:, 0], states[:, 1], pressure_high, pressure_low,
    ))
    shaft_rpm = states[:, 4] * 60 / (2 * np.pi)
    load_resistance = pto.load_controller.resistance(shaft_rpm, states[:, 6])
    body_position = np.column_stack((
        lever * np.sin(angle), np.zeros_like(angle),
        cg_z + lever * (np.cos(angle) - 1),
    ))
    body_velocity = velocity_history[:, :3]
    return OSWECHydraulicResponse(
        time=time, pitch=angle, pitch_velocity=speed,
        body_position=body_position, body_velocity=body_velocity,
        crank_stroke=crank_stroke, crank_speed=crank_jacobian * speed,
        crank_torque=torque, cylinder_force=cylinder_force,
        cylinder_pressure=states[:, :2],
        accumulator_pressure=np.column_stack((pressure_high, pressure_low)),
        valve_flow=valve_flow,
        motor_torque=pto.motor.torque(pressure_high - pressure_low),
        motor_flow=pto.motor.flow(states[:, 4]),
        shaft_speed_rpm=shaft_rpm, generator_current=states[:, 5],
        generator_voltage=load_resistance * states[:, 5],
        load_resistance=load_resistance,
        wave_elevation=sea.elevation, excitation_force=excitation_force,
    )

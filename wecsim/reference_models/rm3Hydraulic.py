"""Coupled two-heave RM3 runner for the published rectified hydraulic PTO."""

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from wecsim.hydrodynamics.bodyClass import BodyClass
from wecsim.components.hydraulic import RectifiedHydraulicPTO
from wecsim.environment.irregularWave import IrregularComponents, synthesize_irregular_response


@dataclass(frozen=True)
class RM3HydraulicResponse:
    time: np.ndarray
    body_heave: np.ndarray  # time, float/spar, world z position
    body_heave_velocity: np.ndarray
    pto_stroke: np.ndarray
    pto_velocity: np.ndarray
    pto_force: np.ndarray
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
    excitation_heave: np.ndarray  # time, float/spar


def run_rm3_rectified_hydraulic(
    hydro_file: str | Path,
    components: IrregularComponents,
    pto: RectifiedHydraulicPTO,
    *,
    dt: float = .01,
    end_time: float = 400.,
    ramp_time: float = 100.,
    radiation_memory: float = 60.,
    initial_displacement: tuple[float, float] = (0., 0.),
    initial_velocity: tuple[float, float] = (0., 0.),
    rho: float = 1000.,
    g: float = 9.81,
) -> RM3HydraulicResponse:
    """Integrate float and spar heave with a discrete hydraulic PTO network.

    This is the vertical-slider layout in ``RM3_cHydraulic_PTO``: both bodies
    are constrained to heave, radiation uses the HDF5 impulse-response
    kernel, and PTO state updates use forward Euler. The body coordinates
    use an implicit trapezoidal step. The provided sea components determine
    the phase realization; use a saved MATLAB phase for an exact input pair.
    """
    values = np.asarray((dt, end_time, ramp_time, radiation_memory, rho, g),
                        dtype=float)
    if (not isinstance(components, IrregularComponents)
            or not isinstance(pto, RectifiedHydraulicPTO)):
        raise TypeError("RM3 hydraulic runner needs a sea and rectified PTO")
    if (not np.isfinite(values).all() or dt <= 0 or end_time <= 0
            or ramp_time < 0 or radiation_memory <= 0 or rho <= 0 or g <= 0):
        raise ValueError("RM3 hydraulic time and fluid settings are invalid")
    steps = round(end_time / dt)
    memory_steps = round(radiation_memory / dt)
    if (not np.isclose(steps * dt, end_time, rtol=0, atol=1e-10)
            or not np.isclose(memory_steps * dt, radiation_memory,
                              rtol=0, atol=1e-10)):
        raise ValueError("duration and memory must be integer time steps")
    initial_q = np.asarray(initial_displacement, dtype=float)
    initial_v = np.asarray(initial_velocity, dtype=float)
    if (initial_q.shape != (2,) or initial_v.shape != (2,)
            or not np.isfinite(initial_q).all()
            or not np.isfinite(initial_v).all()):
        raise ValueError("initial float/spar heaves and speeds must be finite pairs")

    time = np.arange(steps + 1) * dt
    convolution_time = np.arange(memory_steps + 1) * dt
    sea = tuple(synthesize_irregular_response(
        hydro_file, components, dt=dt, end_time=end_time,
        ramp_time=ramp_time, body_number=number, rho=rho, g=g,
    ) for number in (1, 2))
    excitation = np.column_stack([record.excitation_force[:, 2]
                                  for record in sea])

    effective_mass = np.empty(2)
    restoring = np.empty(2)
    static_force = np.empty(2)
    center_heave = np.empty(2)
    kernel = np.empty((2, len(convolution_time)))
    for index in range(2):
        body = BodyClass(str(hydro_file))
        body.bodyNumber = index + 1
        body.bodyTotal = 2
        body.readH5file()
        if int(np.asarray(body.dof).item()) != 6:
            raise ValueError("RM3 hydraulic runner needs two six-DOF BEM bodies")
        irf_time = body.hydroData["hydro_coeffs"]["radiation_damping"][
            "impulse_response_fun"]["t"]
        if radiation_memory > np.max(irf_time) + 1e-10:
            raise ValueError("radiation memory exceeds the HDF5 kernel")
        body.mass = "equilibrium"
        body.hydroStiffness = np.zeros((6, 6))
        body.viscDrag = {"Drag": np.zeros((6, 6)), "cd": np.zeros(6),
                         "characteristicArea": np.zeros(6)}
        body.linearDamping = np.zeros((6, 6))
        body.hydroForcePre(
            [], [0], len(convolution_time), convolution_time, [], dt,
            rho, g, "noWaveCIC", np.vstack((time, np.zeros_like(time))),
            index + 1, 2, 0, 0, 0,
        )
        hydro = body.hydroForce
        rigid_mass = float(np.asarray(body.mass).item())
        effective_mass[index] = rigid_mass + hydro["fAddedMass"][2, 2]
        restoring[index] = hydro["linearHydroRestCoef"][2, 2]
        static_force[index] = (
            rho * float(np.asarray(body.dispVol).item()) - rigid_mass
        ) * g
        center_heave[index] = np.asarray(body.cg).ravel()[2]
        kernel[index] = np.asarray(hydro["irkb"])[:, 2, 2]
    if (not np.isfinite(np.r_[effective_mass, restoring, static_force,
                               center_heave, kernel.ravel()]).all()
            or np.any(effective_mass <= 0)):
        raise ValueError("RM3 heave hydrodynamics are nonfinite or singular")

    q = np.zeros((steps + 1, 2))
    v = np.zeros_like(q)
    a = np.zeros_like(q)
    q[0], v[0] = initial_q, initial_v
    states = np.empty((steps + 1, 7))
    states[0] = pto.initial_state()
    force = np.empty(steps + 1)
    force[0] = pto.force(states[0])
    force_axis = np.array([1., -1.])
    a[0] = (excitation[0] + static_force - restoring * q[0]
            - dt / 2 * kernel[:, 0] * v[0]
            + force[0] * force_axis) / effective_mass

    for step in range(1, steps + 1):
        previous = step - 1
        rates = pto.state_rate(q[previous, 0] - q[previous, 1],
                               v[previous, 0] - v[previous, 1],
                               states[previous])
        states[step] = states[previous] + dt * rates
        if not np.isfinite(states[step]).all():
            raise FloatingPointError("hydraulic PTO state became nonfinite")
        force[step] = pto.force(states[step])

        memory = min(step, memory_steps)
        known_radiation = np.empty(2)
        for index in range(2):
            known_radiation[index] = dt * np.dot(
                kernel[index, 1:memory + 1],
                v[step - memory:step, index][::-1],
            )
            if step >= memory_steps:
                known_radiation[index] -= (
                    dt / 2 * kernel[index, -1] * v[step - memory, index]
                )
        rhs = (excitation[step] + static_force
               - restoring * (q[previous] + dt / 2 * v[previous])
               - known_radiation + force[step] * force_axis)
        v[step] = (
            effective_mass * (v[previous] + dt / 2 * a[previous])
            + dt / 2 * rhs
        ) / (effective_mass + dt**2 / 4 * (restoring + kernel[:, 0]))
        q[step] = q[previous] + dt / 2 * (v[previous] + v[step])
        a[step] = (
            excitation[step] + static_force - restoring * q[step]
            - known_radiation - dt / 2 * kernel[:, 0] * v[step]
            + force[step] * force_axis
        ) / effective_mass

    pressure_high = pto.high_accumulator.pressure(states[:, 2])
    pressure_low = pto.low_accumulator.pressure(states[:, 3])
    valve_flow = np.column_stack(pto.valve.flows(
        states[:, 0], states[:, 1], pressure_high, pressure_low,
    ))
    shaft_rpm = states[:, 4] * 60 / (2 * np.pi)
    load_resistance = pto.load_controller.resistance(shaft_rpm, states[:, 6])
    return RM3HydraulicResponse(
        time=time, body_heave=q + center_heave,
        body_heave_velocity=v, pto_stroke=q[:, 0] - q[:, 1],
        pto_velocity=v[:, 0] - v[:, 1], pto_force=force,
        cylinder_pressure=states[:, :2],
        accumulator_pressure=np.column_stack((pressure_high, pressure_low)),
        valve_flow=valve_flow,
        motor_torque=pto.motor.torque(pressure_high - pressure_low),
        motor_flow=pto.motor.flow(states[:, 4]),
        shaft_speed_rpm=shaft_rpm, generator_current=states[:, 5],
        generator_voltage=load_resistance * states[:, 5],
        load_resistance=load_resistance,
        wave_elevation=sea[0].elevation, excitation_heave=excitation,
    )

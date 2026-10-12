"""Coupled seven-coordinate dynamics for the floating OWC application.

The floater translates and rotates in six world coordinates; the water column
slides along its local z axis. Hydrodynamics use the two uncoupled 6x6 blocks
of a regular-wave BEM file, as in the published WEC-Sim application. MoorDyn,
chamber pressure, and turbine speed advance with the body state at each step.
"""

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.spatial.transform import Rotation

from wecsim.hydrodynamics.bodyClass import BodyClass
from wecsim.components.floatingOwc import (FloatingOwcChamber, FloatingOwcColumnJoint,
                          FloatingOwcTurbine)
from wecsim.dynamics.moorDyn import MoorDyn


@dataclass(frozen=True)
class FloatingOwcResponse:
    """World body states and air-train outputs sampled at every time step.

    Poses and velocities have six columns: world xyz then xyz Euler angles or
    world angular velocity. ``mooring_pose`` and ``mooring_velocity`` are the
    endpoint predictions passed to native MoorDyn at each coupling step.
    ``pto_mechanical_power`` is positive for power entering the PTO and
    includes spring exchange; ``pto_dissipated_power`` is damper loss only.
    """

    time: np.ndarray
    floater_pose: np.ndarray
    floater_velocity: np.ndarray
    column_pose: np.ndarray
    column_velocity: np.ndarray
    stroke: np.ndarray
    stroke_speed: np.ndarray
    pto_force: np.ndarray
    pto_mechanical_power: np.ndarray
    pto_dissipated_power: np.ndarray
    chamber_pressure: np.ndarray
    turbine_speed: np.ndarray
    mooring_pose: np.ndarray
    mooring_velocity: np.ndarray
    mooring_force: np.ndarray
    turbine_power: np.ndarray
    pneumatic_power: np.ndarray


def _cross_matrix(vector):
    x, y, z = vector
    return np.array([[0, -z, y], [z, 0, -x], [-y, x, 0]])


def _body_data(h5_file, number, mass, inertia, frequency, rho, g, dt):
    body = BodyClass(str(h5_file))
    body.bodyNumber = number
    body.bodyTotal = 2
    body.readH5file()
    body.mass = mass
    body.inertia = np.asarray(inertia, dtype=float)
    body.hydroStiffness = np.zeros((6, 6))
    body.viscDrag = {"Drag": np.zeros((6, 6)), "cd": np.zeros(6),
                     "characteristicArea": np.zeros(6)}
    body.linearDamping = np.zeros((6, 6))
    body.hydroForcePre(
        frequency, [0], 1, np.array([0.0]), [], dt, rho, g,
        "regular", np.zeros((2, 2)), number, 2, 0, 0, 0,
    )
    hydro = body.hydroForce
    center = np.asarray(body.cg, dtype=float).reshape(3)
    buoyancy_center = np.asarray(body.cb, dtype=float).reshape(3)
    volume = float(np.asarray(body.dispVol).item())
    mass_value = float(np.asarray(body.mass).item())
    rigid_mass = np.diag([mass_value] * 3 + list(inertia))
    restoring_at_equilibrium = np.r_[
        np.zeros(3), -np.cross(buoyancy_center - center,
                               [0, 0, rho * g * volume]),
    ]
    return {
        "center": center,
        "mass": rigid_mass + np.asarray(hydro["fAddedMass"]),
        "damping": np.asarray(hydro["fDamping"]),
        "restoring": np.asarray(hydro["linearHydroRestCoef"]),
        "restoring_at_equilibrium": restoring_at_equilibrium,
        "vertical_bias": (rho * volume - mass_value) * g,
        "excitation_real": np.asarray(hydro["fExt"]["re"]),
        "excitation_imag": np.asarray(hydro["fExt"]["im"]),
    }


def solve_floating_owc(
    h5_file: str | Path,
    moordyn: MoorDyn,
    *,
    dt: float = 0.01,
    end_time: float = 500.0,
    wave_height: float = 4.5,
    wave_period: float = 11.2,
    ramp_time: float = 50.0,
    floater_mass: str | float = "equilibrium",
    floater_inertia: tuple[float, float, float] = (1.531e9, 1.531e9, 0.1118e9),
    column_mass: float = 4_493_450.0,
    column_height: float = 50.69,
    column_diameter: float = 5.89,
    column_inertia: tuple[float, float, float] | None = None,
    pto_stiffness: float = 0.0,
    pto_damping: float = 0.0,
    moordyn_point: tuple[float, float, float] | None = None,
    chamber: FloatingOwcChamber | None = None,
    turbine: FloatingOwcTurbine | None = None,
    initial_state: np.ndarray | None = None,
    initial_rotor_speed: float = 150.0,
    rho: float = 1000.0,
    g: float = 9.81,
) -> FloatingOwcResponse:
    """Advance body, native mooring, chamber, and turbine together.

    Defaults reproduce the published ``OWC/FloatingOWC`` regular-wave case.
    ``moordyn_point`` is the body-local coupling point relative to floater CG;
    by default it is the world origin in equilibrium. ``initial_state``, when
    provided, is floater xyz, xyz Euler angles, slider stroke, floater world
    linear/angular velocity, slider speed, gauge pressure, and rotor speed.
    Positive PTO stiffness and damping oppose the slider stroke and speed;
    the returned ``pto_force`` acts on the water column along the slider.
    The native MoorDyn model is stepped once per fixed interval using a
    predicted endpoint. Other state derivatives use RK4 on that interval.
    """
    if not isinstance(moordyn, MoorDyn):
        raise TypeError("moordyn must be a MoorDyn session")
    scalars = (dt, end_time, wave_period, wave_height, ramp_time,
               column_mass, column_height, column_diameter,
               pto_stiffness, pto_damping, rho, g,
               initial_rotor_speed)
    if (not np.isfinite(scalars).all() or dt <= 0 or end_time <= 0
            or wave_period <= 0 or wave_height < 0 or ramp_time < 0
            or min(column_mass, column_height, column_diameter, rho, g,
                   initial_rotor_speed) <= 0):
        raise ValueError("floating OWC physical parameters must be finite and valid")
    steps = round(end_time / dt)
    if not np.isclose(steps * dt, end_time, rtol=0, atol=1e-9):
        raise ValueError("end_time must be an integer number of time steps")
    floater_inertia = np.asarray(floater_inertia, dtype=float)
    if column_inertia is None:
        radius = column_diameter / 2
        transverse = column_mass * (3 * radius**2 + column_height**2) / 12
        column_inertia = (transverse, transverse, column_mass * radius**2 / 2)
    column_inertia = np.asarray(column_inertia, dtype=float)
    if (floater_inertia.shape != (3,) or column_inertia.shape != (3,)
            or not np.isfinite(floater_inertia).all()
            or not np.isfinite(column_inertia).all()
            or np.any(floater_inertia <= 0) or np.any(column_inertia <= 0)):
        raise ValueError("body inertia must be three finite positive values")

    frequency = 2 * np.pi / wave_period
    bodies = (
        _body_data(h5_file, 1, floater_mass, floater_inertia,
                   frequency, rho, g, dt),
        _body_data(h5_file, 2, column_mass, column_inertia,
                   frequency, rho, g, dt),
    )
    separation = bodies[1]["center"][2] - bodies[0]["center"][2]
    joint = FloatingOwcColumnJoint(separation)
    point = np.asarray(-bodies[0]["center"] if moordyn_point is None
                       else moordyn_point, dtype=float)
    if point.shape != (3,) or not np.isfinite(point).all():
        raise ValueError("moordyn_point must be a finite body-local three-vector")
    if chamber is None:
        area = np.pi * column_diameter**2 / 4
        chamber = FloatingOwcChamber(area, area * 4.5, 1.4,
                                     101325, 1.25, 0.75, 0.775)
    if turbine is None:
        turbine = FloatingOwcTurbine()
    if initial_state is None:
        state = np.r_[bodies[0]["center"], np.zeros(12), initial_rotor_speed]
    else:
        state = np.asarray(initial_state, dtype=float).copy()
    if state.shape != (16,) or not np.isfinite(state).all() or state[15] <= 0:
        raise ValueError("initial_state must be 16 finite values with positive rotor speed")

    floater_jacobian = np.zeros((6, 7))
    floater_jacobian[:, :6] = np.eye(6)
    time = np.arange(steps + 1) * dt
    history = np.empty((steps + 1, 16))
    history[0] = state
    mooring_pose = np.empty((steps + 1, 6))
    mooring_velocity = np.empty((steps + 1, 6))
    mooring_force = np.zeros((steps + 1, 6))

    def connection(values):
        rotated = Rotation.from_euler("xyz", values[3:6]).apply(point)
        pose = np.r_[values[:3] + rotated, values[3:6]]
        velocity = np.r_[values[7:10] + np.cross(values[10:13], rotated),
                         values[10:13]]
        jacobian = np.zeros((6, 7))
        jacobian[:3, :3] = np.eye(3)
        jacobian[:3, 3:6] = -_cross_matrix(rotated)
        jacobian[3:6, 3:6] = np.eye(3)
        return pose, velocity, jacobian

    def derivative(at_time, values, line_force):
        position = values[:3]
        angles = values[3:6]
        stroke = values[6]
        speed = values[7:14]
        axis = Rotation.from_euler("xyz", angles).apply([0, 0, 1])
        lever = (separation + stroke) * axis
        column_jacobian = np.zeros((6, 7))
        column_jacobian[:3, :3] = np.eye(3)
        column_jacobian[:3, 3:6] = -_cross_matrix(lever)
        column_jacobian[:3, 6] = axis
        column_jacobian[3:6, 3:6] = np.eye(3)
        angular = speed[3:6]
        axis_rate = np.cross(angular, axis)
        column_bias = np.r_[
            2 * speed[6] * axis_rate
            + (separation + stroke) * np.cross(angular, axis_rate),
            np.zeros(3),
        ]
        mass = np.zeros((7, 7))
        force = np.zeros(7)
        ramp = (1 if ramp_time == 0 or at_time >= ramp_time
                else (1 - np.cos(np.pi * at_time / ramp_time)) / 2)
        for body, jacobian, bias, center_position in (
            (bodies[0], floater_jacobian, np.zeros(6), position),
            (bodies[1], column_jacobian, column_bias, position + lever),
        ):
            displacement = np.r_[center_position - body["center"], angles]
            velocity = jacobian @ speed
            excitation = wave_height / 2 * ramp * (
                body["excitation_real"] * np.cos(frequency * at_time)
                - body["excitation_imag"] * np.sin(frequency * at_time)
            )
            body_force = (excitation - body["damping"] @ velocity
                          - body["restoring"] @ displacement
                          - body["restoring_at_equilibrium"]
                          - body["mass"] @ bias)
            body_force[2] += body["vertical_bias"]
            mass += jacobian.T @ body["mass"] @ jacobian
            force += jacobian.T @ body_force
        _, _, mooring_jacobian = connection(values)
        force += mooring_jacobian.T @ line_force
        force[6] += float(chamber.force_on_column(values[14]))
        force[6] -= pto_stiffness * stroke + pto_damping * speed[6]
        acceleration = np.linalg.solve(mass, force)

        pitch, yaw = angles[1], angles[2]
        euler_to_angular = np.array([
            [np.cos(yaw) * np.cos(pitch), -np.sin(yaw), 0],
            [np.sin(yaw) * np.cos(pitch), np.cos(yaw), 0],
            [-np.sin(pitch), 0, 1],
        ])
        angle_rate = np.linalg.solve(euler_to_angular, angular)
        column_velocity = column_jacobian @ speed
        column_heave = position[2] + lever[2] - bodies[1]["center"][2]
        pressure_rate = chamber.pressure_derivative(
            values[14], column_heave, column_velocity[2], values[15])
        turbine_rate = turbine.evaluate(values[14], values[15]).speed_derivative
        return np.r_[speed[:3], angle_rate, speed[6], acceleration,
                     pressure_rate, turbine_rate]

    mooring_pose[0], mooring_velocity[0], _ = connection(state)
    with moordyn.start(mooring_pose[0], mooring_velocity[0]):
        for index in range(steps):
            at_time = time[index]
            prior_force = mooring_force[index]
            prediction = state + dt * derivative(at_time, state, prior_force)
            pose, velocity, _ = connection(prediction)
            mooring_pose[index + 1] = pose
            mooring_velocity[index + 1] = velocity
            line_force = moordyn.step(pose, velocity, at_time, dt)
            mooring_force[index + 1] = line_force
            k1 = derivative(at_time, state, line_force)
            k2 = derivative(at_time + dt / 2, state + dt * k1 / 2, line_force)
            k3 = derivative(at_time + dt / 2, state + dt * k2 / 2, line_force)
            k4 = derivative(at_time + dt, state + dt * k3, line_force)
            state = state + dt / 6 * (k1 + 2*k2 + 2*k3 + k4)
            if not np.isfinite(state).all():
                raise ValueError("floating OWC integration became nonfinite")
            history[index + 1] = state

    floater_pose = history[:, :6]
    floater_velocity = history[:, 7:13]
    column_pose = joint.column_pose(floater_pose, history[:, 6])
    column_velocity = joint.column_velocity(
        floater_pose, floater_velocity, history[:, 6], history[:, 13])
    power = turbine.evaluate(history[:, 14], history[:, 15])
    pto_force = -pto_stiffness * history[:, 6] - pto_damping * history[:, 13]
    return FloatingOwcResponse(
        time, floater_pose, floater_velocity, column_pose, column_velocity,
        history[:, 6], history[:, 13], pto_force,
        -pto_force * history[:, 13], pto_damping * history[:, 13]**2,
        history[:, 14], history[:, 15],
        mooring_pose, mooring_velocity, mooring_force,
        power.load_power, power.pneumatic_power,
    )

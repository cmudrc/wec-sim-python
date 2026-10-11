"""Reduced four-coordinate model of the published RM3 regular-wave case.

The two bodies share surge and pitch at the floating joint and have separate
heave coordinates. Hydrodynamic coefficients come from production ``BodyClass``
preprocessing. This is a reference-case model, not the general WEC-Sim runner.
"""

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .bodyClass import BodyClass
from .generalDynamics import BodyMotion, DynamicBody, GeneralizedDynamics
from .hardStops import LinearHardStops
from .moorDyn import MoorDyn
from .radiation import RadiationStateSpace


@dataclass(frozen=True)
class RM3RegularResponse:
    """Motion and PTO signals; mechanical power includes spring exchange.

    Positive mechanical power enters the PTO. Dissipated power counts the
    ordinary linear PTO damper and excludes any stop damper.
    """
    time: np.ndarray
    body_position: np.ndarray
    body_velocity: np.ndarray
    coordinate_position: np.ndarray
    coordinate_velocity: np.ndarray
    pto_force: np.ndarray
    pto_stroke: np.ndarray
    pto_velocity: np.ndarray
    pto_mechanical_power: np.ndarray
    pto_dissipated_power: np.ndarray
    pto_stop_force: np.ndarray | None = None
    mooring_surge_position: np.ndarray | None = None
    mooring_surge_force: np.ndarray | None = None
    moordyn_connection_position: np.ndarray | None = None
    moordyn_connection_velocity: np.ndarray | None = None
    moordyn_connection_force: np.ndarray | None = None


def _rm3_slider_motion(coordinate, speed, *, body_index, lever):
    """Map the pitched floating joint and one body-local slider to six DOFs."""
    angle = coordinate[3]
    slide = coordinate[body_index + 1]
    radius = lever + slide
    sine, cosine = np.sin(angle), np.cos(angle)
    jacobian = np.zeros((6, 4))
    jacobian[0, 0] = 1.0
    jacobian[0, body_index + 1] = sine
    jacobian[2, body_index + 1] = cosine
    jacobian[0, 3] = radius * cosine
    jacobian[2, 3] = -radius * sine
    jacobian[4, 3] = 1.0
    curvature = np.array([
        2 * cosine * speed[body_index + 1] * speed[3]
        - radius * sine * speed[3]**2,
        0.0,
        -2 * sine * speed[body_index + 1] * speed[3]
        - radius * cosine * speed[3]**2,
        0.0, 0.0, 0.0,
    ])
    displacement = np.array([
        coordinate[0] + radius * sine, 0.0,
        slide * cosine + lever * (cosine - 1.0),
        0.0, angle, 0.0,
    ])
    return BodyMotion(displacement, jacobian, curvature)


def solve_rm3_regular(
    h5_file: str | Path,
    *,
    wave_height: float = 2.5,
    wave_period: float = 8.0,
    pitch_inertias: tuple[float, float] = (21_306_090.66, 94_407_091.24),
    pto_damping: float = 1_200_000.0,
    pto_stiffness: float = 0.0,
    pto_equilibrium: float = 0.0,
    pto_hard_stops: LinearHardStops | None = None,
    mooring_surge_stiffness: float = 0.0,
    moordyn: MoorDyn | None = None,
    moordyn_point: tuple[float, float, float] | None = None,
    b2b: bool = False,
    radiation_memory: float | None = None,
    radiation_method: str | None = None,
    added_mass_scheme: str = "implicit",
    excitation_force: np.ndarray | None = None,
    no_wave: bool = False,
    initial_coordinate: np.ndarray | None = None,
    initial_speed: np.ndarray | None = None,
    joint_z: float = 0.0,
    dt: float = 0.1,
    end_time: float = 400.0,
    ramp_time: float = 100.0,
    rho: float = 1000.0,
    g: float = 9.81,
) -> RM3RegularResponse:
    """Solve RM3 surge, two heaves, and shared pitch with a relative heave PTO.

    The bodies use equilibrium displaced-volume masses, fixed-frequency
    added mass and radiation damping, hydrostatic restoring, and regular-wave
    excitation. Rigid-body rotation and slider travel change their
    surge/heave Jacobians at each step. Without ``radiation_memory``,
    fixed-frequency damping uses
    classical RK4. With radiation memory, infinite-frequency added mass and
    the radiation impulse-response kernel use a trapezoidal history step.
    ``added_mass_scheme="implicit"`` includes infinite-frequency added mass
    in the effective inertia. ``"simulink_delay"`` reproduces the pinned
    WEC-Sim model's rigid-body mass split and 1e-7 s delayed acceleration
    feedback for numerical comparison. The delay is not a WEC property.
    ``radiation_method="fir"`` instead samples the same kernel as a discrete
    FIR filter and holds its force through each RK4 step.
    ``radiation_method="state_space"`` advances the pinned HDF5 fitted
    radiation states with the body in RK4. The source sets direct term D to
    zero. This fit is nonpassive in active RM3 modes and is opt-in for
    numerical source comparison; convolution remains the physical default.
    A supplied ``pto_hard_stops`` adds unilateral stroke limits to the linear
    PTO and selects adaptive integration with implicit added mass. This path
    currently supports constant-frequency radiation only.
    ``excitation_force`` supplies a sampled two-body, six-DOF wave force at
    every output time for imported irregular spectra; RK4 stage forces are
    linearly interpolated between those samples.
    ``mooring_surge_stiffness`` connects a linear surge spring to the floating
    joint, as in the published RM3 MooringMatrix case.
    ``moordyn`` instead couples a native MoorDyn model at ``moordyn_point``
    in spar-local coordinates. MoorDyn receives a predicted pose and velocity
    once per step; returned loads act during the body solve. The recorded
    connection pose is the pose supplied to MoorDyn.
    ``no_wave=True`` uses the noWaveCIC preprocessing and requires a radiation
    memory. Initial coordinates are shared surge, float heave, spar heave,
    and shared pitch. ``b2b=True`` includes cross-body radiation blocks.
    Sway, roll, yaw, and
    full Simscape joint forces are outside this reduced model.
    """
    inputs = [wave_height, wave_period, pto_damping, pto_stiffness,
              pto_equilibrium, mooring_surge_stiffness,
              joint_z, dt, end_time, ramp_time, rho, g, *pitch_inertias]
    if radiation_memory is not None:
        inputs.append(radiation_memory)
    if not np.isfinite(inputs).all():
        raise ValueError("solver inputs must be finite")
    if (len(pitch_inertias) != 2 or wave_height < 0 or wave_period <= 0
            or any(inertia <= 0 for inertia in pitch_inertias)
            or pto_damping < 0 or pto_stiffness < 0
            or mooring_surge_stiffness < 0 or dt <= 0
            or end_time < 0 or ramp_time < 0 or rho <= 0 or g <= 0):
        raise ValueError("invalid RM3 wave, body, PTO, or time parameters")
    if radiation_memory is not None and radiation_memory <= 0:
        raise ValueError("radiation_memory must be positive")
    if radiation_method is None:
        radiation_method = ("convolution" if radiation_memory is not None
                            else "constant")
    if radiation_method not in ("constant", "convolution", "fir", "state_space"):
        raise ValueError("unsupported RM3 radiation method")
    if radiation_method == "state_space" and radiation_memory is None:
        # BodyClass still needs a small CIC grid for A_inf and excitation;
        # the fitted state evolution itself does not use that memory length.
        radiation_memory = dt
    if added_mass_scheme not in ("implicit", "simulink_delay"):
        raise ValueError("added_mass_scheme must be implicit or simulink_delay")
    if added_mass_scheme == "simulink_delay" and radiation_method != "convolution":
        raise ValueError("simulink_delay requires convolution radiation")
    if ((radiation_method == "constant" and radiation_memory is not None)
            or (radiation_method != "constant" and radiation_memory is None)):
        raise ValueError("constant radiation has no memory; convolution and FIR need it")
    if not isinstance(no_wave, bool):
        raise ValueError("no_wave must be a boolean")
    if no_wave and (wave_height != 0 or radiation_memory is None):
        raise ValueError("no_wave needs zero wave height and radiation memory")
    if not isinstance(b2b, bool):
        raise ValueError("b2b must be a boolean")
    if pto_hard_stops is not None:
        if not isinstance(pto_hard_stops, LinearHardStops):
            raise TypeError("pto_hard_stops must be LinearHardStops")
        if radiation_method != "constant" or added_mass_scheme != "implicit":
            raise ValueError("PTO hard stops currently need constant radiation and implicit added mass")
    if moordyn is not None:
        if not isinstance(moordyn, MoorDyn) or moordyn_point is None:
            raise ValueError("MoorDyn coupling needs a MoorDyn session and spar-local point")
        point = np.asarray(moordyn_point, dtype=float)
        if point.shape != (3,) or not np.isfinite(point).all():
            raise ValueError("moordyn_point must be a finite spar-local three-vector")
        if (radiation_method != "convolution" or pto_hard_stops is not None
                or mooring_surge_stiffness):
            raise ValueError("MoorDyn coupling needs convolution radiation and no other mooring or hard stops")
    elif moordyn_point is not None:
        raise ValueError("moordyn_point needs a MoorDyn session")
    steps = round(end_time / dt)
    if not np.isclose(steps * dt, end_time, rtol=0, atol=1e-10):
        raise ValueError("end_time must be an integer multiple of dt")
    time = np.arange(steps + 1) * dt
    if excitation_force is not None:
        excitation_force = np.asarray(excitation_force, dtype=float)
        if (excitation_force.shape != (steps + 1, 2, 6)
                or not np.isfinite(excitation_force).all()):
            raise ValueError("excitation_force must have shape (time, 2, 6) and be finite")
        if no_wave or radiation_memory is None or wave_height != 0:
            raise ValueError("sampled excitation needs zero regular-wave height and radiation memory")
    omega = 2 * np.pi / wave_period
    if radiation_memory is not None:
        memory_steps = round(radiation_memory / dt)
        if not np.isclose(memory_steps * dt, radiation_memory, rtol=0, atol=1e-10):
            raise ValueError("radiation_memory must be an integer multiple of dt")
        convolution_time = np.arange(memory_steps + 1) * dt
    data = []
    for index, pitch_inertia in enumerate(pitch_inertias, start=1):
        body = BodyClass(str(h5_file))
        body.bodyNumber = index
        body.bodyTotal = 2
        body.readH5file()
        if int(np.asarray(body.dof).item()) != 6:
            raise ValueError("each RM3 hydrodynamic body must have six DOFs")
        body.mass = "equilibrium"
        body.hydroStiffness = np.zeros((6, 6))
        body.viscDrag = {
            "Drag": np.zeros((6, 6)),
            "cd": np.zeros(6),
            "characteristicArea": np.zeros(6),
        }
        body.linearDamping = np.zeros((6, 6))
        if radiation_memory is None:
            body.hydroForcePre(
                omega, [0], 1, np.array([0.0]), [], dt, rho, g,
                "regular", np.vstack((time, np.zeros_like(time))),
                index, 2, 0, 0, int(b2b),
            )
        else:
            irf_time = body.hydroData["hydro_coeffs"]["radiation_damping"][
                "impulse_response_fun"]["t"]
            if radiation_memory > np.max(irf_time) + 1e-10:
                raise ValueError("radiation_memory exceeds the HDF5 kernel")
            body.hydroForcePre(
                [] if no_wave or excitation_force is not None else omega,
                [0], len(convolution_time),
                convolution_time, [], dt, rho, g,
                "noWaveCIC" if no_wave or excitation_force is not None
                else "regularCIC",
                np.vstack((time, np.zeros_like(time))),
                index, 2, 0, 0, int(b2b),
            )
        mass = float(np.asarray(body.mass).item())
        center = np.asarray(body.cg).ravel()
        if not np.isclose(center[0:2], 0.0, atol=1e-10).all():
            raise ValueError("the RM3 model requires body centers on the joint axis")
        lever = float(center[2] - joint_z)
        rigid_mass = np.diag([mass, mass, mass, 0.0, pitch_inertia, 0.0])
        hydro = body.hydroForce
        if b2b:
            added_mass = tuple(np.asarray(hydro["fAddedMass"])[:, 6*j:6*(j+1)] for j in range(2))
            damping = (tuple(np.asarray(hydro["fDamping"])[:, 6*j:6*(j+1)] for j in range(2))
                       if radiation_memory is None else tuple(np.zeros((6, 6)) for _ in range(2)))
        else:
            added_mass = [np.zeros((6, 6)), np.zeros((6, 6))]
            damping = [np.zeros((6, 6)), np.zeros((6, 6))]
            added_mass[index - 1] = np.asarray(hydro["fAddedMass"])
            if radiation_memory is None:
                damping[index - 1] = np.asarray(hydro["fDamping"])
        if radiation_method in ("constant", "state_space"):
            kernel = None
        else:
            raw_kernel = np.asarray(hydro["irkb"])
            if b2b:
                kernel = raw_kernel
            else:
                kernel = np.zeros((len(raw_kernel), 6, 12))
                kernel[:, :, 6 * (index - 1):6 * index] = raw_kernel
        data.append({
            "center_z": float(center[2]),
            "lever": lever,
            "rigid_mass": rigid_mass,
            "added_mass": added_mass,
            "damping": damping,
            "kernel": kernel,
            "restoring": np.asarray(hydro["linearHydroRestCoef"]),
            "re": np.asarray(hydro["fExt"]["re"]),
            "im": np.asarray(hydro["fExt"]["im"]),
            "vertical_bias": (rho * float(np.asarray(body.dispVol).item()) - mass) * g,
        })

    pto_coupling = np.zeros((4, 4))
    pto_coupling[1, 1] = pto_coupling[2, 2] = 1
    pto_coupling[1, 2] = pto_coupling[2, 1] = -1
    mechanical_stiffness = pto_stiffness * pto_coupling
    mechanical_stiffness[0, 0] = mooring_surge_stiffness

    dynamic_bodies = []
    for index, body in enumerate(data):
        lever = body["lever"]

        def motion(q, v, *, index=index, lever=lever):
            return _rm3_slider_motion(q, v, body_index=index, lever=lever)

        def excitation(at_time, *, re=body["re"], im=body["im"],
                       body_index=index):
            if excitation_force is not None:
                sample = min(max(at_time / dt, 0.0), float(steps))
                left = min(int(sample), steps - 1)
                fraction = sample - left
                return ((1 - fraction) * excitation_force[left, body_index]
                        + fraction * excitation_force[left + 1, body_index])
            ramp = (1.0 if ramp_time == 0 or at_time >= ramp_time
                    else (1.0 - np.cos(np.pi * at_time / ramp_time)) / 2)
            return wave_height / 2 * ramp * (
                re * np.cos(omega * at_time) - im * np.sin(omega * at_time)
            )

        dynamic_bodies.append(DynamicBody(
            rigid_mass=body["rigid_mass"],
            added_mass=tuple(body["added_mass"]),
            damping=tuple(body["damping"]),
            restoring=body["restoring"],
            static_force=np.array([0, 0, body["vertical_bias"], 0, 0, 0]),
            reference_position=np.array([0, 0, body["center_z"], 0, 0, 0]),
            motion=motion,
            excitation=excitation,
            radiation_kernel=body["kernel"],
        ))

    def stop_generalized_force(coordinate, speed):
        stroke = coordinate[1] - coordinate[2]
        stroke_speed = speed[1] - speed[2]
        reaction = float(pto_hard_stops.force(stroke, stroke_speed))
        return np.array([0.0, reaction, -reaction, 0.0])

    system = GeneralizedDynamics(
        tuple(dynamic_bodies), 4,
        pto_stiffness=mechanical_stiffness,
        pto_damping=pto_damping * pto_coupling,
        pto_equilibrium=np.array([0, pto_equilibrium / 2,
                                  -pto_equilibrium / 2, 0]),
        radiation_discretization=("fir" if radiation_method == "fir"
                                  else "trapezoid"),
        added_mass_delay=(1e-7 if added_mass_scheme == "simulink_delay" else None),
        nonlinear_force=(stop_generalized_force
                         if pto_hard_stops is not None else None),
    )
    fitted_radiation = (
        RadiationStateSpace.from_rm3_hdf5(h5_file, body_to_body=b2b)
        if radiation_method == "state_space" else None
    )
    connection_position = connection_velocity = connection_force = None
    if moordyn is not None:
        connection_position = np.zeros((steps + 1, 6))
        connection_velocity = np.zeros_like(connection_position)
        connection_force = np.zeros_like(connection_position)
        spar = dynamic_bodies[1]

        def connection_motion(q, v):
            body_motion = spar.motion(q, v)
            angle = q[3]
            sine, cosine = np.sin(angle), np.cos(angle)
            rotated = np.array([
                cosine * point[0] + sine * point[2], point[1],
                -sine * point[0] + cosine * point[2],
            ])
            tangent = np.array([
                -sine * point[0] + cosine * point[2], 0.0,
                -cosine * point[0] - sine * point[2],
            ])
            jacobian = body_motion.jacobian.copy()
            jacobian[:3, 3] += tangent
            position = spar.reference_position + body_motion.displacement
            position[:3] += rotated
            return position, jacobian @ v, jacobian

        q0 = (np.zeros(4) if initial_coordinate is None
              else np.asarray(initial_coordinate, dtype=float))
        v0 = (np.zeros(4) if initial_speed is None
              else np.asarray(initial_speed, dtype=float))
        connection_position[0], connection_velocity[0], _ = connection_motion(q0, v0)

        def step_force(at_time, interval, predicted_coordinate, predicted_speed):
            index = round(at_time / dt) + 1
            position, velocity, jacobian = connection_motion(
                predicted_coordinate, predicted_speed)
            connection_position[index] = position
            connection_velocity[index] = velocity
            connection_force[index] = moordyn.step(
                position, velocity, at_time, interval)
            return jacobian.T @ connection_force[index]

        with moordyn.start(connection_position[0], connection_velocity[0]):
            solved = system.integrate(
                dt=dt, end_time=end_time,
                initial_coordinate=initial_coordinate, initial_speed=initial_speed,
                step_force=step_force,
            )
    else:
        solved = system.integrate(
            dt=dt, end_time=end_time,
            initial_coordinate=initial_coordinate, initial_speed=initial_speed,
            adaptive_regular=pto_hard_stops is not None,
            radiation_state_space=fitted_radiation,
        )
    q = solved.coordinate
    v = solved.speed
    pto_stroke = q[:, 1] - q[:, 2]
    pto_velocity = v[:, 1] - v[:, 2]
    pto_stop_force = (pto_hard_stops.force(pto_stroke, pto_velocity)
                      if pto_hard_stops is not None else None)
    pto_force = (-pto_damping * pto_velocity
                 - pto_stiffness * (pto_stroke - pto_equilibrium)
                 + (pto_stop_force if pto_stop_force is not None else 0))
    return RM3RegularResponse(
        time=solved.time, body_position=solved.body_position,
        body_velocity=solved.body_velocity,
        coordinate_position=q, coordinate_velocity=v,
        pto_force=pto_force, pto_stroke=pto_stroke,
        pto_velocity=pto_velocity,
        pto_mechanical_power=-pto_force * pto_velocity,
        pto_dissipated_power=pto_damping * pto_velocity**2,
        pto_stop_force=pto_stop_force,
        mooring_surge_position=q[:, 0] if mooring_surge_stiffness else None,
        mooring_surge_force=(-mooring_surge_stiffness * q[:, 0]
                             if mooring_surge_stiffness else None),
        moordyn_connection_position=connection_position,
        moordyn_connection_velocity=connection_velocity,
        moordyn_connection_force=connection_force,
    )

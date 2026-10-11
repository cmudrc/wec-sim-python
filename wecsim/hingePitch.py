"""Forced pitch response of one hydrodynamic body on a horizontal hinge.

The solver accepts a six-component excitation force history in the body's
equilibrium global axes. Its hydrodynamic coefficients come from the production
``BodyClass`` preprocessing. A separate wave model must supply excitation;
this module does not synthesize irregular waves.
"""

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import trimesh

from .bodyClass import BodyClass
from .generalDynamics import BodyMotion, DynamicBody, GeneralizedDynamics
from .morison import finite_depth_wavenumber
from .nonlinearHydro import mesh_pressure_wrenches, regular_wave_mesh_pressures


@dataclass(frozen=True)
class HingePitchResponse:
    time: np.ndarray
    angle: np.ndarray
    angular_velocity: np.ndarray
    center_position: np.ndarray
    center_velocity: np.ndarray
    excitation_torque: np.ndarray
    pto_torque: np.ndarray
    excitation_force: np.ndarray | None = None


def solve_hinged_pitch_from_excitation(
    h5_file: str | Path,
    excitation_force: np.ndarray,
    *,
    hinge_z: float,
    body_mass: float,
    pitch_inertia: float,
    pto_damping: float,
    pto_stiffness: float = 0.0,
    pto_equilibrium: float = 0.0,
    dt: float = 0.1,
    memory_time: float = 30.0,
    added_mass_scheme: str = "implicit",
    rho: float = 1000.0,
    g: float = 9.81,
) -> HingePitchResponse:
    """Integrate hinge pitch with hydrodynamic radiation and a linear PTO.

    ``excitation_force`` is an ``(N, 6)`` array with columns surge, sway,
    heave, roll, pitch, yaw at uniform ``dt``. The body starts at its HDF5
    equilibrium position, zero angle, and zero velocity. The hinge is parallel
    to the global y axis. The model includes nonlinear rigid-body kinematics,
    hydrostatic restoring, infinite-frequency added mass, and the radiation
    impulse-response convolution. It assumes a fixed hinge and no body-to-body
    hydrodynamic interaction, Morison forces, or other DOFs.

    The fixed-step trapezoidal update uses a finite memory window. It matches
    the published OSWEC example's 0.1 s step and 30 s radiation memory when
    supplied with that run's excitation forces.
    """
    force = np.asarray(excitation_force, dtype=float)
    if force.ndim != 2 or force.shape[1] != 6 or len(force) < 2:
        raise ValueError("excitation_force must have shape (N, 6), N >= 2")
    if not np.isfinite(force).all():
        raise ValueError("excitation_force must be finite")
    values = [hinge_z, body_mass, pitch_inertia, pto_damping,
              pto_stiffness, pto_equilibrium, dt, memory_time, rho, g]
    if not np.isfinite(values).all():
        raise ValueError("solver parameters must be finite")
    if (body_mass <= 0 or pitch_inertia <= 0 or pto_damping < 0
            or pto_stiffness < 0 or dt <= 0 or memory_time <= 0
            or rho <= 0 or g <= 0):
        raise ValueError("mass, inertia, dt, memory time, rho, and g must be positive")
    if added_mass_scheme not in ("implicit", "simulink_delay"):
        raise ValueError("added_mass_scheme must be implicit or simulink_delay")
    memory_steps = round(memory_time / dt)
    if not np.isclose(memory_steps * dt, memory_time, rtol=0, atol=1e-10):
        raise ValueError("memory_time must be an integer multiple of dt")

    body = BodyClass(str(h5_file))
    body.bodyNumber = 1
    body.readH5file()
    if int(np.asarray(body.dof).item()) != 6:
        raise ValueError("the hinged-pitch solver requires one six-DOF body")
    cg_z = float(body.hydroData["properties"]["cg"][0, 2])
    lever = cg_z - hinge_z
    if lever <= 0:
        raise ValueError("the hinge must be below the center of gravity")
    irf_time = body.hydroData["hydro_coeffs"]["radiation_damping"]["impulse_response_fun"]["t"]
    if memory_time > np.max(irf_time) + 1e-10:
        raise ValueError("memory_time exceeds the radiation kernel in the HDF5 file")

    count = len(force)
    time = np.arange(count) * dt
    convolution_time = np.arange(memory_steps + 1) * dt
    body.mass = body_mass
    body.hydroStiffness = np.zeros((6, 6))
    body.viscDrag = {
        "Drag": np.zeros((6, 6)),
        "cd": np.zeros(6),
        "characteristicArea": np.zeros(6),
    }
    body.linearDamping = np.zeros((6, 6))
    body.hydroForcePre(
        [], [0], len(convolution_time), convolution_time, [], dt, rho, g,
        "noWaveCIC", np.vstack((time, np.zeros_like(time))),
        1, 1, 0, 0, 0,
    )

    displaced_volume = float(np.asarray(body.dispVol).item())
    static_force = np.array([
        0.0, 0.0, (rho * displaced_volume - body_mass) * g, 0.0, 0.0, 0.0,
    ])

    def motion(q, v):
        angle, speed = q[0], v[0]
        j = np.array([
            lever * np.cos(angle), 0.0, -lever * np.sin(angle),
            0.0, 1.0, 0.0,
        ])[:, None]
        bias = np.array([
            -lever * np.sin(angle), 0.0, -lever * np.cos(angle),
            0.0, 0.0, 0.0,
        ]) * speed**2
        displacement = np.array([
            lever * np.sin(angle), 0.0,
            lever * (np.cos(angle) - 1.0), 0.0, angle, 0.0,
        ])
        return BodyMotion(displacement, j, bias)

    def excitation(at_time):
        index = round(at_time / dt)
        if index < 0 or index >= count or not np.isclose(index * dt, at_time, atol=1e-9):
            raise ValueError("sampled excitation is defined only on the simulation time grid")
        return force[index]

    rigid_mass = np.diag([
        body_mass, body_mass, body_mass, 0.0, pitch_inertia, 0.0,
    ])
    device = GeneralizedDynamics((DynamicBody(
        rigid_mass=rigid_mass,
        added_mass=(np.asarray(body.hydroForce["fAddedMass"]),),
        damping=(np.zeros((6, 6)),),
        restoring=np.asarray(body.hydroForce["linearHydroRestCoef"]),
        static_force=static_force,
        reference_position=np.array([0, 0, cg_z, 0, 0, 0]),
        motion=motion,
        excitation=excitation,
        radiation_kernel=np.asarray(body.hydroForce["irkb"]),
    ),), 1,
        pto_stiffness=np.array([[pto_stiffness]]),
        pto_damping=np.array([[pto_damping]]),
        pto_equilibrium=np.array([pto_equilibrium]),
        added_mass_delay=(1e-7 if added_mass_scheme == "simulink_delay" else None),
    )
    solved = device.integrate(dt=dt, end_time=(count - 1) * dt)
    angle = solved.coordinate[:, 0]
    speed = solved.speed[:, 0]
    center_position = solved.body_position[:, 0, :3]
    center_velocity = solved.body_velocity[:, 0, :3]
    excitation_torque = np.array([
        motion(solved.coordinate[i], solved.speed[i]).jacobian[:, 0] @ force[i]
        for i in range(count)
    ])
    return HingePitchResponse(
        time=time,
        angle=angle,
        angular_velocity=speed,
        center_position=center_position,
        center_velocity=center_velocity,
        excitation_torque=excitation_torque,
        pto_torque=(-pto_damping * speed
                    - pto_stiffness * (angle - pto_equilibrium)),
        excitation_force=force,
    )


def solve_hinged_pitch_regular(
    h5_file: str | Path,
    *,
    wave_height: float,
    wave_period: float,
    hinge_z: float,
    body_mass: float,
    pitch_inertia: float,
    pto_damping: float,
    pto_stiffness: float = 0.0,
    pto_equilibrium: float = 0.0,
    dt: float = 0.1,
    end_time: float = 400.0,
    ramp_time: float = 100.0,
    rho: float = 1000.0,
    g: float = 9.81,
) -> HingePitchResponse:
    """Integrate the hinged flap under a regular wave and frequency damping."""
    values = [wave_height, wave_period, hinge_z, body_mass, pitch_inertia,
              pto_damping, pto_stiffness, pto_equilibrium, dt, end_time,
              ramp_time, rho, g]
    if not np.isfinite(values).all():
        raise ValueError("solver parameters must be finite")
    if (wave_height < 0 or wave_period <= 0 or body_mass <= 0
            or pitch_inertia <= 0 or pto_damping < 0 or pto_stiffness < 0
            or dt <= 0 or end_time < 0 or ramp_time < 0 or rho <= 0 or g <= 0):
        raise ValueError("invalid wave, body, PTO, or time parameters")
    steps = round(end_time / dt)
    if not np.isclose(steps * dt, end_time, rtol=0, atol=1e-10):
        raise ValueError("end_time must be an integer multiple of dt")
    time = np.arange(steps + 1) * dt
    omega = 2 * np.pi / wave_period

    body = BodyClass(str(h5_file))
    body.bodyNumber = 1
    body.bodyTotal = 1
    body.readH5file()
    if int(np.asarray(body.dof).item()) != 6:
        raise ValueError("the hinged flap needs six hydrodynamic DOFs")
    cg_z = float(body.hydroData["properties"]["cg"][0, 2])
    lever = cg_z - hinge_z
    if lever <= 0:
        raise ValueError("the hinge must be below the center of gravity")
    body.mass = body_mass
    body.hydroStiffness = np.zeros((6, 6))
    body.viscDrag = {
        "Drag": np.zeros((6, 6)), "cd": np.zeros(6),
        "characteristicArea": np.zeros(6),
    }
    body.linearDamping = np.zeros((6, 6))
    body.hydroForcePre(
        omega, [0], 1, np.array([0.0]), [], dt, rho, g,
        "regular", np.vstack((time, np.zeros_like(time))),
        1, 1, 0, 0, 0,
    )
    hydro = body.hydroForce
    re = np.asarray(hydro["fExt"]["re"])
    im = np.asarray(hydro["fExt"]["im"])

    def motion(q, v):
        angle, speed = q[0], v[0]
        jacobian = np.array([
            lever * np.cos(angle), 0.0, -lever * np.sin(angle),
            0.0, 1.0, 0.0,
        ])[:, None]
        bias = np.array([
            -lever * np.sin(angle), 0.0, -lever * np.cos(angle),
            0.0, 0.0, 0.0,
        ]) * speed**2
        displacement = np.array([
            lever * np.sin(angle), 0.0,
            lever * (np.cos(angle) - 1.0), 0.0, angle, 0.0,
        ])
        return BodyMotion(displacement, jacobian, bias)

    def excitation(at_time):
        ramp = (1.0 if ramp_time == 0 or at_time >= ramp_time
                else (1.0 - np.cos(np.pi * at_time / ramp_time)) / 2)
        return wave_height / 2 * ramp * (
            re * np.cos(omega * at_time) - im * np.sin(omega * at_time)
        )

    rigid_mass = np.diag([
        body_mass, body_mass, body_mass, 0.0, pitch_inertia, 0.0,
    ])
    vertical_bias = (
        rho * float(np.asarray(body.dispVol).item()) - body_mass
    ) * g
    device = GeneralizedDynamics((DynamicBody(
        rigid_mass=rigid_mass,
        added_mass=(np.asarray(hydro["fAddedMass"]),),
        damping=(np.asarray(hydro["fDamping"]),),
        restoring=np.asarray(hydro["linearHydroRestCoef"]),
        static_force=np.array([0, 0, vertical_bias, 0, 0, 0]),
        reference_position=np.array([0, 0, cg_z, 0, 0, 0]),
        motion=motion, excitation=excitation,
    ),), 1,
        pto_stiffness=np.array([[pto_stiffness]]),
        pto_damping=np.array([[pto_damping]]),
        pto_equilibrium=np.array([pto_equilibrium]),
    )
    solved = device.integrate(dt=dt, end_time=end_time)
    angle = solved.coordinate[:, 0]
    speed = solved.speed[:, 0]
    force_history = np.array([excitation(at_time) for at_time in time])
    torque = np.array([
        motion(solved.coordinate[i], solved.speed[i]).jacobian[:, 0]
        @ force_history[i] for i in range(len(time))
    ])
    return HingePitchResponse(
        time=time, angle=angle, angular_velocity=speed,
        center_position=solved.body_position[:, 0, :3],
        center_velocity=solved.body_velocity[:, 0, :3],
        excitation_torque=torque,
        pto_torque=(-pto_damping * speed
                    - pto_stiffness * (angle - pto_equilibrium)),
        excitation_force=force_history,
    )


def solve_hinged_pitch_nonlinear_regular(
    h5_file: str | Path,
    geometry_file: str | Path,
    *,
    wave_height: float,
    wave_period: float,
    hinge_z: float,
    body_mass: float,
    pitch_inertia: float,
    pto_damping: float,
    pto_stiffness: float = 0.0,
    pto_equilibrium: float = 0.0,
    dt: float = 0.1,
    end_time: float = 120.0,
    ramp_time: float = 40.0,
    rho: float = 1000.0,
    g: float = 9.81,
) -> HingePitchResponse:
    """Advance a hinged flap with instantaneous regular-wave mesh pressures.

    The pressure and frequency-domain force laws are paired with the pinned
    nonlinear OSWEC case. The body motion uses implicit added mass and resolves
    pressure forces between output samples; it does not adopt Simulink's delayed
    acceleration feedback.
    """
    values = [wave_height, wave_period, hinge_z, body_mass, pitch_inertia,
              pto_damping, pto_stiffness, pto_equilibrium, dt, end_time,
              ramp_time, rho, g]
    if (not np.isfinite(values).all() or wave_height < 0 or wave_period <= 0
            or body_mass <= 0 or pitch_inertia <= 0 or pto_damping < 0
            or pto_stiffness < 0 or dt <= 0 or end_time < 0
            or ramp_time < 0 or rho <= 0 or g <= 0):
        raise ValueError("invalid nonlinear hinged-flap parameters")
    steps = round(end_time / dt)
    if not np.isclose(steps * dt, end_time, rtol=0, atol=1e-10):
        raise ValueError("end_time must be an integer multiple of dt")
    time = np.arange(steps + 1) * dt
    omega = 2 * np.pi / wave_period

    body = BodyClass(str(h5_file))
    body.bodyNumber = body.bodyTotal = 1
    body.readH5file()
    if int(np.asarray(body.dof).item()) != 6:
        raise ValueError("the hinged flap needs six hydrodynamic DOFs")
    cg = np.asarray(body.hydroData["properties"]["cg"], dtype=float).ravel()
    if not np.allclose(cg[:2], 0, atol=1e-10, rtol=0):
        raise ValueError("nonlinear hinged flap needs its CG on the hinge's z axis")
    lever = cg[2] - hinge_z
    if lever <= 0:
        raise ValueError("the hinge must be below the center of gravity")
    depth = float(np.asarray(
        body.hydroData["simulation_parameters"]["water_depth"]
    ).item())
    if not np.isfinite(depth) or depth <= 0:
        raise ValueError("nonlinear regular waves need finite positive water depth")
    wave_number = finite_depth_wavenumber(
        np.array([omega]), water_depth=depth, gravity=g,
    )[0]
    mesh = trimesh.load_mesh(str(geometry_file), process=False)
    if len(mesh.faces) == 0:
        raise ValueError("nonlinear hinged flap needs a triangular geometry mesh")

    body.mass = body_mass
    body.hydroStiffness = np.zeros((6, 6))
    body.viscDrag = {
        "Drag": np.zeros((6, 6)), "cd": np.zeros(6),
        "characteristicArea": np.zeros(6),
    }
    body.linearDamping = np.zeros((6, 6))
    body.hydroForcePre(
        omega, [0], 1, np.array([0.0]), [], dt, rho, g,
        "regular", np.vstack((time, np.zeros_like(time))),
        1, 1, 0, 2, 0,
    )
    hydro = body.hydroForce
    re = np.asarray(hydro["fExt"]["re"])
    im = np.asarray(hydro["fExt"]["im"])

    def motion(q, v):
        angle, speed = q[0], v[0]
        jacobian = np.array([
            lever * np.cos(angle), 0.0, -lever * np.sin(angle),
            0.0, 1.0, 0.0,
        ])[:, None]
        bias = np.array([
            -lever * np.sin(angle), 0.0, -lever * np.cos(angle),
            0.0, 0.0, 0.0,
        ]) * speed**2
        displacement = np.array([
            lever * np.sin(angle), 0.0,
            lever * (np.cos(angle) - 1), 0.0, angle, 0.0,
        ])
        return BodyMotion(displacement, jacobian, bias)

    def forces(at_time, angle):
        pose = np.array([[
            lever * np.sin(angle), 0.0,
            hinge_z + lever * np.cos(angle), 0.0, angle, 0.0,
        ]])
        pressure = regular_wave_mesh_pressures(
            mesh.vertices, mesh.faces, pose, [at_time],
            center_gravity=cg, rho=rho, gravity=g,
            water_depth=depth, wave_number=wave_number,
            wave_height=wave_height, wave_period=wave_period,
            ramp_time=ramp_time,
        )
        wrench = mesh_pressure_wrenches(
            mesh.vertices, mesh.faces, pose, pressure,
        )
        ramp = (1.0 if ramp_time == 0 or at_time >= ramp_time
                else (1 - np.cos(np.pi * at_time / ramp_time)) / 2)
        excitation = (wave_height / 2 * ramp * (
            re * np.cos(omega * at_time) - im * np.sin(omega * at_time)
        ) + ramp * (wrench.nonlinear_wave[0] - wrench.linear_wave[0]))
        hydrostatic = wrench.hydrostatic[0].copy()
        hydrostatic[2] -= body_mass * g
        return hydrostatic, excitation

    def state_excitation(at_time, coordinate, speed):
        hydrostatic, excitation = forces(at_time, coordinate[0])
        return hydrostatic + excitation

    rigid_mass = np.diag([
        body_mass, body_mass, body_mass, 0.0, pitch_inertia, 0.0,
    ])
    device = GeneralizedDynamics((DynamicBody(
        rigid_mass=rigid_mass,
        added_mass=(np.asarray(hydro["fAddedMass"]),),
        damping=(np.asarray(hydro["fDamping"]),),
        restoring=np.zeros((6, 6)), static_force=np.zeros(6),
        reference_position=np.r_[cg, np.zeros(3)],
        motion=motion, excitation=lambda at_time: np.zeros(6),
        state_excitation=state_excitation,
    ),), 1,
        pto_stiffness=np.array([[pto_stiffness]]),
        pto_damping=np.array([[pto_damping]]),
        pto_equilibrium=np.array([pto_equilibrium]),
    )
    solved = device.integrate(
        dt=dt, end_time=end_time, adaptive_regular=True,
    )
    angle = solved.coordinate[:, 0]
    speed = solved.speed[:, 0]
    excitation_force = np.array([
        forces(at_time, state_angle)[1]
        for at_time, state_angle in zip(solved.time, angle)
    ])
    torque = np.array([
        motion(solved.coordinate[i], solved.speed[i]).jacobian[:, 0]
        @ excitation_force[i]
        for i in range(len(solved.time))
    ])
    return HingePitchResponse(
        time=solved.time, angle=angle, angular_velocity=speed,
        center_position=solved.body_position[:, 0, :3],
        center_velocity=solved.body_velocity[:, 0, :3],
        excitation_torque=torque,
        pto_torque=(-pto_damping * speed
                    - pto_stiffness * (angle - pto_equilibrium)),
        excitation_force=excitation_force,
    )

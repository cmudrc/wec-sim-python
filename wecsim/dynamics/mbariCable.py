"""Coupled planar dynamics of the published MBARI Cable application."""

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from wecsim.hydrodynamics.bodyClass import BodyClass
from wecsim.dynamics.cable import PlanarCableAttachment, WecSimCableTension
from wecsim.dynamics.generalDynamics import BodyMotion, DynamicBody, GeneralizedDynamics


@dataclass(frozen=True)
class MBARICableResponse:
    time: np.ndarray
    wave_elevation: np.ndarray
    body_position: np.ndarray
    body_velocity: np.ndarray
    cable_displacement: np.ndarray
    cable_speed: np.ndarray
    cable_force_z: np.ndarray


def _body_motion(index, q, v):
    jacobian = np.zeros((6, 7))
    bias = np.zeros(6)
    if index < 2:
        start = 0 if index == 0 else 3
        jacobian[0, start] = jacobian[2, start + 1] = 1
        jacobian[4, start + 2] = 1
        displacement = np.zeros(6)
        displacement[[0, 2, 4]] = q[start:start + 3]
    else:
        buoy_pitch, cylinder_pitch = q[2], q[6]
        buoy_rate, cylinder_rate = v[2], v[6]
        jacobian[0, 0] = jacobian[2, 1] = jacobian[4, 6] = 1
        jacobian[0, 2] = -.69 * np.cos(buoy_pitch)
        jacobian[0, 6] = -4.21 * np.cos(cylinder_pitch)
        jacobian[2, 2] = .69 * np.sin(buoy_pitch)
        jacobian[2, 6] = 4.21 * np.sin(cylinder_pitch)
        displacement = np.array([
            q[0] - .69 * np.sin(buoy_pitch) - 4.21 * np.sin(cylinder_pitch),
            0.,
            q[1] + .69 * (1 - np.cos(buoy_pitch))
            + 4.21 * (1 - np.cos(cylinder_pitch)),
            0., cylinder_pitch, 0.,
        ])
        bias[0] = (.69 * np.sin(buoy_pitch) * buoy_rate**2
                   + 4.21 * np.sin(cylinder_pitch) * cylinder_rate**2)
        bias[2] = (.69 * np.cos(buoy_pitch) * buoy_rate**2
                   + 4.21 * np.cos(cylinder_pitch) * cylinder_rate**2)
    return BodyMotion(displacement, jacobian, bias)


def run_mbari_cable(
    hydro_file: str | Path,
    cable: WecSimCableTension | None = None,
    attachment: PlanarCableAttachment | None = None,
    *,
    wave_height: float = 1.,
    wave_period: float = 8.,
    dt: float = .01,
    end_time: float = 120.,
    ramp_time: float = 40.,
    rho: float = 1000.,
    g: float = 9.81,
    endpoint_mass: float = 1.,
) -> MBARICableResponse:
    """Advance all three bodies with cable tension and endpoint drag.

    The buoy and cylinder share the published spherical joint at z=-0.59 m.
    The spar and cylinder are connected by a spring cable with two drag bodies
    at its endpoints. Each published drag body contributes 1 kg of
    translational inertia at its body-local attachment point. The endpoint
    mass can be changed explicitly while the cable and drag force laws remain
    the same.
    Hydrodynamics use the published regular-wave HDF5 bodies and implicit
    frequency-dependent added mass. This is a planar case runner, not a
    general six-DOF Simscape constraint solver.
    """
    if cable is None:
        cable = WecSimCableTension(1_000_000., 100., 17.8, 18.)
    if attachment is None:
        attachment = PlanarCableAttachment((0., 1.95), (0., -5.2), 18.)
    if not isinstance(cable, WecSimCableTension) or not isinstance(
            attachment, PlanarCableAttachment):
        raise TypeError("cable and attachment need their public configuration objects")
    values = np.asarray((wave_height, wave_period, dt, end_time,
                         ramp_time, rho, g, endpoint_mass), dtype=float)
    if (not np.isfinite(values).all() or wave_height < 0 or wave_period <= 0
            or dt <= 0 or end_time <= 0 or ramp_time < 0 or rho <= 0 or g <= 0
            or endpoint_mass < 0
            or not np.isclose(cable.initial_length, attachment.initial_length,
                              rtol=0, atol=1e-10)
            or attachment.base_offset[0] != 0
            or attachment.follower_offset[0] != 0):
        raise ValueError("MBARI wave, cable, and simulation settings are invalid")
    steps = round(end_time / dt)
    if not np.isclose(steps * dt, end_time, rtol=0, atol=1e-10):
        raise ValueError("duration must have an integer number of steps")

    time = np.arange(steps + 1) * dt
    omega = 2 * np.pi / wave_period
    centers = (
        np.array([0., 0., .1, 0., 0., 0.]),
        np.array([0., 0., -29.95, 0., 0., 0.]),
        np.array([0., 0., -4.8, 0., 0., 0.]),
    )
    masses = (1080., 815., 600.)
    pitch_inertias = (541.3, 815 * (3 * 1.4**2 + .1**2) / 12, 3898.4)
    drag_cd = (
        (1.15, 1.15, 1., .5, .5, 0.),
        (.8, .8, 0., 1., 1., 0.),
        (1.15, 1.15, 0., 1.15, 1.15, 0.),
    )
    drag_area = (
        (2.9568, 2.9568, 5.4739, 5.4739, 5.4739, 0.),
        (1.12, 1.12, 6.56, 6.56, 6.56, 0.),
        (1.496, 1.496, 0., 1.496, 1.496, 0.),
    )

    def motion(index):
        return lambda q, v: _body_motion(index, q, v)

    def ramp(at_time):
        return (1. if ramp_time == 0 or at_time >= ramp_time else
                (1 - np.cos(np.pi * at_time / ramp_time)) / 2)

    def endpoint_inertia(body_index, offset):
        """Return a point mass's wrench and spatial mass at a rotating point."""
        angle_index = 5 if body_index == 1 else 6
        local_x, local_z = offset

        def evaluate(_at_time, q, v):
            angle = q[angle_index]
            rate = v[angle_index]
            sine, cosine = np.sin(angle), np.cos(angle)
            arm_x = local_x * cosine + local_z * sine
            arm_z = -local_x * sine + local_z * cosine
            point_jacobian = np.zeros((2, 6))
            point_jacobian[0, 0] = point_jacobian[1, 2] = 1
            point_jacobian[0, 4] = arm_z
            point_jacobian[1, 4] = -arm_x
            relative_bias = -np.array([arm_x, arm_z]) * rate**2
            mass = endpoint_mass * point_jacobian.T @ point_jacobian
            # Source joint loads contain endpoint m*a but no separate weight.
            force = -endpoint_mass * point_jacobian.T @ relative_bias
            return force, mass

        return evaluate

    zero = np.zeros((6, 6))
    bodies = []
    for index in range(3):
        if index < 2:
            hydro_body = BodyClass(str(hydro_file))
            hydro_body.bodyNumber = index + 1
            hydro_body.bodyTotal = 2
            hydro_body.readH5file()
            if (int(np.asarray(hydro_body.dof).item()) != 6
                    or not np.allclose(np.asarray(hydro_body.cg).ravel(),
                                       centers[index][:3], rtol=0, atol=1e-8)):
                raise ValueError("MBARI HDF5 body geometry does not match the published case")
            hydro_body.mass = masses[index]
            hydro_body.hydroStiffness = zero.copy()
            hydro_body.viscDrag = {
                "Drag": zero.copy(), "cd": np.zeros(6),
                "characteristicArea": np.zeros(6),
            }
            hydro_body.linearDamping = zero.copy()
            hydro_body.hydroForcePre(
                omega, [0], 1, np.array([0.]), [], dt, rho, g,
                "regular", np.vstack((time, np.zeros_like(time))),
                index + 1, 2, 0, 0, 0,
            )
            hydro = hydro_body.hydroForce
            added = np.asarray(hydro["fAddedMass"])
            damping = np.asarray(hydro["fDamping"])
            restoring = np.asarray(hydro["linearHydroRestCoef"])
            re = np.asarray(hydro["fExt"]["re"])
            im = np.asarray(hydro["fExt"]["im"])
            volume = float(np.asarray(hydro_body.dispVol).item())
        else:
            added = damping = restoring = zero
            re = im = np.zeros(6)
            volume = .2
        added_blocks = [zero.copy() for _ in range(3)]
        damping_blocks = [zero.copy() for _ in range(3)]
        added_blocks[index] = added
        damping_blocks[index] = damping
        rigid = np.diag([masses[index]] * 3
                        + [0., pitch_inertias[index], 0.])

        def make_excitation(body_index, real, imaginary):
            coefficients = .5 * wave_height
            body_drag = (.5 * rho * np.asarray(drag_cd[body_index])
                         * np.asarray(drag_area[body_index]))

            def excitation(at_time):
                return coefficients * ramp(at_time) * (
                    real * np.cos(omega * at_time)
                    - imaginary * np.sin(omega * at_time)
                )

            def state_force(at_time, q, v):
                velocity = _body_motion(body_index, q, v).jacobian @ v
                return excitation(at_time) - body_drag * np.abs(velocity) * velocity

            return excitation, state_force

        excitation, state_force = make_excitation(index, re, im)
        point_inertia = None
        if endpoint_mass and index in (1, 2):
            offset = (attachment.base_offset if index == 1
                      else attachment.follower_offset)
            point_inertia = endpoint_inertia(index, offset)
        bodies.append(DynamicBody(
            rigid_mass=rigid,
            added_mass=tuple(added_blocks),
            damping=tuple(damping_blocks),
            restoring=restoring,
            static_force=np.array([0., 0., (rho * volume - masses[index]) * g,
                                   0., 0., 0.]),
            reference_position=centers[index],
            motion=motion(index),
            excitation=excitation,
            state_excitation=state_force,
            state_inertia=point_inertia,
        ))

    def generalized_cable_force(q, v):
        base_motion = _body_motion(1, q, v)
        follower_motion = _body_motion(2, q, v)
        base_pose = centers[1] + base_motion.displacement
        follower_pose = centers[2] + follower_motion.displacement
        base_rate = base_motion.jacobian @ v
        follower_rate = follower_motion.jacobian @ v
        axes = (0, 2, 4)
        displacement, speed = attachment.motion(
            base_pose[list(axes)], base_rate[list(axes)],
            follower_pose[list(axes)], follower_rate[list(axes)],
        )
        force = cable.force_z(displacement, speed)
        base_wrench, follower_wrench = attachment.world_wrenches(
            base_pose[list(axes)], follower_pose[list(axes)], force,
        )

        def add_endpoint_drag(wrench, pose, rate, offset):
            local_x, local_z = offset
            sine, cosine = np.sin(pose[4]), np.cos(pose[4])
            arm_x = local_x * cosine + local_z * sine
            arm_z = -local_x * sine + local_z * cosine
            velocity = np.array([rate[0] + arm_z * rate[4],
                                 rate[2] - arm_x * rate[4]])
            drag = -.5 * rho * 1.4 * 10 * np.abs(velocity) * velocity
            wrench = wrench.copy()
            wrench[:2] += drag
            wrench[2] += arm_z * drag[0] - arm_x * drag[1]
            return np.array([wrench[0], 0., wrench[1], 0., wrench[2], 0.])

        base_load = add_endpoint_drag(
            base_wrench, base_pose, base_rate, attachment.base_offset,
        )
        follower_load = add_endpoint_drag(
            follower_wrench, follower_pose, follower_rate,
            attachment.follower_offset,
        )
        return (base_motion.jacobian.T @ base_load
                + follower_motion.jacobian.T @ follower_load)

    solved = GeneralizedDynamics(
        tuple(bodies), 7, nonlinear_force=generalized_cable_force,
    ).integrate(dt=dt, end_time=end_time, adaptive_regular=True)
    base_pose = solved.body_position[:, 1][:, (0, 2, 4)]
    base_rate = solved.body_velocity[:, 1][:, (0, 2, 4)]
    follower_pose = solved.body_position[:, 2][:, (0, 2, 4)]
    follower_rate = solved.body_velocity[:, 2][:, (0, 2, 4)]
    displacement, speed = attachment.motion(
        base_pose, base_rate, follower_pose, follower_rate,
    )
    wave_ramp = np.array([ramp(at_time) for at_time in time])
    return MBARICableResponse(
        time=time,
        wave_elevation=.5 * wave_height * wave_ramp * np.cos(omega * time),
        body_position=solved.body_position,
        body_velocity=solved.body_velocity,
        cable_displacement=displacement,
        cable_speed=speed,
        cable_force_z=cable.force_z(displacement, speed),
    )

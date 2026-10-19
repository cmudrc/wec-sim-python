"""Opt-in published-step integration for one OSWEC pitch hinge.

WEC-Sim breaks its added-mass feedback loop with a very short Transport Delay.
This path approximates its stage-wise delayed signal from the two preceding
accepted body accelerations. The resulting trajectory is step-sensitive; the
physical default continues to solve added mass implicitly.
"""

from dataclasses import replace

import numpy as np

from wecsim.dynamics.generalDynamics import DynamicsResponse, GeneralizedDynamics


def integrate_nonlinear_hinge_ode4(device: GeneralizedDynamics, *,
                                   dt: float, end_time: float,
                                   delay: float = 1e-7,
                                   initial_angle: float = 0.0,
                                   initial_angular_velocity: float = 0.0) -> DynamicsResponse:
    """Advance one regular-wave hinge with source-style delayed inertia.

    The excitation may depend only on time or on the current flap pose.
    """
    if (len(device.bodies) != 1 or device.coordinate_count != 1
            or device.bodies[0].radiation_kernel is not None
            or device.bodies[0].state_inertia is not None
            or device.nonlinear_force is not None
            or device.controlled_ptos or device.linear_generators
            or device.added_mass_delay is not None):
        raise ValueError("source hinge ode4 needs one continuous regular-wave hinge")
    if (not np.isfinite([dt, end_time, delay, initial_angle,
                         initial_angular_velocity]).all() or dt <= 0
            or end_time < 0 or delay <= 0 or delay >= dt):
        raise ValueError("source hinge ode4 needs finite time and a shorter delay")
    steps = round(end_time / dt)
    if not np.isclose(steps * dt, end_time, rtol=0, atol=1e-10):
        raise ValueError("end_time must be an integer multiple of dt")

    body = device.bodies[0]
    added = body.added_mass[0]
    # WEC-Sim places this rigid mass shift in the inertia matrix and applies
    # only the remaining added mass through its delayed force path.
    shifted = np.zeros((6, 6))
    shifted[:3, :3] = 2 * np.trace(added[:3, :3]) * np.eye(3)
    for axis in range(3, 6):
        shifted[axis, axis] = added[axis, axis]
        for other in range(axis + 1, 6):
            shifted[axis, other] = shifted[other, axis] = added[axis, other]
    applied_added = added - shifted
    source_body = replace(
        body, rigid_mass=body.rigid_mass + shifted,
        added_mass=(np.zeros((6, 6)),),
    )
    source = GeneralizedDynamics(
        (source_body,), 1, pto_stiffness=device.pto_stiffness,
        pto_damping=device.pto_damping,
        pto_equilibrium=device.pto_equilibrium, pto_bias=device.pto_bias,
    )

    time = np.arange(steps + 1) * dt
    coordinate = np.zeros((steps + 1, 1))
    speed = np.zeros_like(coordinate)
    coordinate[0, 0] = initial_angle
    speed[0, 0] = initial_angular_velocity
    acceleration = np.zeros_like(coordinate)
    previous_body_accelerations = []

    def accelerate(at_time, q, v, delayed):
        motion = body.motion(q, v)
        delayed_force = motion.jacobian.T @ (-applied_added @ delayed)
        return source.acceleration(at_time, q, v,
                                   applied_force=delayed_force)

    for step, at_time in enumerate(time):
        if step == 0:
            current_delay = np.zeros(6)
        elif step == 1:
            current_delay = previous_body_accelerations[0]
        else:
            previous, last = previous_body_accelerations[-2:]
            current_delay = last + (1 - delay / dt) * (last - previous)
        q, v = coordinate[step], speed[step]
        acceleration[step] = accelerate(at_time, q, v, current_delay)
        motion = body.motion(q, v)
        previous_body_accelerations.append(
            motion.jacobian @ acceleration[step] + motion.bias_acceleration
        )
        if step == steps:
            break

        # The delay at a minor RK4 stage uses accepted major-step signals;
        # it never feeds that stage's trial acceleration back into itself.
        def derivative(stage_time, state, fraction):
            if fraction == 0:
                delayed = current_delay
            elif step == 0:
                delayed = previous_body_accelerations[-1]
            else:
                previous, last = previous_body_accelerations[-2:]
                delayed = last + (fraction - delay / dt) * (last - previous)
            stage_q, stage_v = state[:1], state[1:]
            return np.r_[stage_v, accelerate(stage_time, stage_q,
                                             stage_v, delayed)]

        state = np.r_[q, v]
        k1 = derivative(at_time, state, 0)
        k2 = derivative(at_time + dt / 2, state + dt * k1 / 2, .5)
        k3 = derivative(at_time + dt / 2, state + dt * k2 / 2, .5)
        k4 = derivative(at_time + dt, state + dt * k3, 1)
        advanced = state + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6
        if not np.isfinite(advanced).all():
            raise ValueError("source hinge ode4 produced nonfinite motion")
        coordinate[step + 1], speed[step + 1] = advanced[:1], advanced[1:]

    position = np.zeros((steps + 1, 1, 6))
    velocity = np.zeros_like(position)
    for step in range(steps + 1):
        motion = body.motion(coordinate[step], speed[step])
        position[step, 0] = body.reference_position + motion.displacement
        velocity[step, 0] = motion.jacobian @ speed[step]
    return DynamicsResponse(time, coordinate, speed, acceleration,
                            position, velocity)

"""Fixed-step WEC-Sim convolution and delayed-mass integration.

The source evaluates its convolution at accepted major steps and extrapolates
three saved radiation forces with PCHIP during ode4's minor stages. Its tiny
added-mass delay likewise uses accepted accelerations, with the value at the
first stage carried forward from the preceding step.
"""

import numpy as np
from scipy.interpolate import PchipInterpolator


def integrate_source_cic_ode4(device, time, coordinate, speed, acceleration,
                              dt):
    """Advance a two-body CIC system through the source's numerical path."""
    bodies = device.bodies
    count = len(bodies)
    kernels = [body.radiation_kernel for body in bodies]
    velocity_history = np.zeros((len(time), count * 6))
    body_acceleration_history = []
    radiation_history = []
    zero_acceleration = np.zeros((count, 6))
    accepted_delay = zero_acceleration
    n = device.coordinate_count

    def save_velocity(index):
        for body_index, body in enumerate(bodies):
            motion = body.motion(coordinate[index], speed[index])
            velocity_history[index, 6 * body_index:6 * (body_index + 1)] = (
                motion.jacobian @ speed[index]
            )

    def radiation(index):
        forces = []
        for kernel in kernels:
            memory = min(index, len(kernel) - 1)
            force = dt * np.einsum(
                "tij,tj->i", kernel[1:memory + 1],
                velocity_history[index - memory:index][::-1],
            )
            if len(kernel) > 1 and index >= len(kernel) - 1:
                force -= dt / 2 * (
                    kernel[-1] @ velocity_history[index - memory]
                )
            force += dt / 2 * (kernel[0] @ velocity_history[index])
            forces.append(force)
        return tuple(forces)

    def accelerate(at_time, q, v, radiation_force, delayed_acceleration):
        device.delayed_body_acceleration = tuple(delayed_acceleration)
        return device.acceleration(
            at_time, q, v, known_radiation=radiation_force, dt=dt,
            complete_radiation=True,
        )

    def save_body_acceleration(index):
        values = []
        for body in bodies:
            motion = body.motion(coordinate[index], speed[index])
            values.append(
                motion.jacobian @ acceleration[index]
                + motion.bias_acceleration
            )
        body_acceleration_history.append(np.stack(values))

    save_velocity(0)
    radiation_history.append(radiation(0))
    acceleration[0] = accelerate(
        time[0], coordinate[0], speed[0], radiation_history[0],
        zero_acceleration,
    )
    save_body_acceleration(0)

    for index in range(len(time) - 1):
        # The source's persistent convolution block starts with two zero
        # samples at negative times. PCHIP then predicts minor-step loads.
        recent = radiation_history[-3:]
        values = np.asarray(recent)
        if len(recent) < 3:
            values = np.concatenate((
                np.zeros((3 - len(recent), count, 6)), values,
            ))
        extrapolator = PchipInterpolator(
            [-2, -1, 0], values, axis=0,
            extrapolate=True,
        )
        radiation_half = tuple(extrapolator(.5))
        radiation_end = tuple(extrapolator(1))

        def delayed(fraction):
            if fraction == 0:
                return accepted_delay
            if index == 0:
                return body_acceleration_history[-1]
            previous, last = body_acceleration_history[-2:]
            return last + (fraction - device.added_mass_delay / dt) * (
                last - previous
            )

        def derivative(at_time, state, fraction, radiation_force):
            q, v = state[:n], state[n:]
            return np.r_[v, accelerate(
                at_time, q, v, radiation_force, delayed(fraction),
            )]

        state = np.r_[coordinate[index], speed[index]]
        at_time = time[index]
        k1 = derivative(at_time, state, 0, radiation_history[-1])
        k2 = derivative(at_time + dt / 2, state + dt * k1 / 2,
                        .5, radiation_half)
        k3 = derivative(at_time + dt / 2, state + dt * k2 / 2,
                        .5, radiation_half)
        k4 = derivative(at_time + dt, state + dt * k3, 1, radiation_end)
        advanced = state + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6
        if not np.isfinite(advanced).all():
            raise ValueError("source CIC ode4 produced nonfinite motion")
        coordinate[index + 1], speed[index + 1] = advanced[:n], advanced[n:]

        save_velocity(index + 1)
        radiation_history.append(radiation(index + 1))
        next_delay = delayed(1)
        acceleration[index + 1] = accelerate(
            time[index + 1], coordinate[index + 1], speed[index + 1],
            radiation_history[-1], next_delay,
        )
        accepted_delay = next_delay
        save_body_acceleration(index + 1)

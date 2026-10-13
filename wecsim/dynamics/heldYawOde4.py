"""Published WEC-Sim ode4 path for a single held passive-yaw coordinate.

The source samples convolution at major steps, extrapolates its last three
loads to RK4 minor stages, and delays its selected excitation heading by one
major step. This opt-in path preserves those numerical details without
changing the default implicit radiation solver.
"""

import numpy as np
from scipy.interpolate import PchipInterpolator

from wecsim.hydrodynamics.passiveYaw import HeldPassiveYawExcitation


def integrate_held_yaw_ode4(system, time, q, v, a, dt):
    """Fill one-coordinate motion arrays using the source's sampled step."""
    moving = [index for index, body in enumerate(system.bodies)
              if isinstance(body.state_excitation, HeldPassiveYawExcitation)]
    if (system.coordinate_count != 1 or len(moving) != 1
            or system.nonlinear_force is not None
            or system.radiation_discretization != "trapezoid"
            or any(body.state_inertia is not None
                   or body.radiation_kernel is None for body in system.bodies)):
        raise ValueError("wecsim_ode4 needs one pure held-yaw coordinate")
    moving_index = moving[0]
    held = system.bodies[moving_index].state_excitation
    for index, body in enumerate(system.bodies):
        motion = body.motion(np.zeros(1), np.zeros(1))
        expected = np.zeros((6, 1))
        if index == moving_index:
            expected[5, 0] = 1
        if (not np.array_equal(motion.jacobian, expected)
                or not np.array_equal(motion.bias_acceleration, np.zeros(6))):
            raise ValueError("wecsim_ode4 needs pure yaw and stationary other bodies")
    if (len(held.model.time) != 2 * len(time) - 1
            or (len(time) > 1 and not np.isclose(
                held.model.time[1], dt / 2, rtol=0, atol=1e-10))):
        raise ValueError("wecsim_ode4 needs half-step sampled excitation")

    # The fixed bodies have zero velocity and zero generalized Jacobian, so
    # only the moving body's yaw-to-yaw IRF enters the scalar equation.
    kernel = system.bodies[moving_index].radiation_kernel[
        :, 5, 6 * moving_index + 5,
    ]
    radiation_times = np.array([-3 * dt, -2 * dt, -dt])
    radiation_values = np.zeros(3)

    def acceleration(at_time, angle, speed, load):
        known = [np.zeros(6) for _ in system.bodies]
        known[moving_index][5] = load
        return system.acceleration(
            at_time, np.array([angle]), np.array([speed]),
            known_radiation=tuple(known), dt=dt,
            complete_radiation=True,
        )[0]

    for step, at_time in enumerate(time):
        memory = min(step, len(kernel) - 1)
        history = v[step - memory:step + 1, 0][::-1]
        radiation = dt * np.dot(kernel[:memory + 1], history)
        radiation -= dt / 2 * kernel[0] * v[step, 0]
        if step >= len(kernel) - 1:
            radiation -= dt / 2 * kernel[-1] * v[step - memory, 0]
        radiation_times = np.r_[radiation_times[1:], at_time]
        radiation_values = np.r_[radiation_values[1:], radiation]

        # The source's selected heading is logged at the current major step,
        # but its Memory block exposes the previous selection to RK4 stages.
        current_heading = held._heading(float(q[step, 0]))
        held.force_history.append(held._force(
            at_time, float(q[step, 0]), current_heading,
        ))
        a[step, 0] = acceleration(at_time, q[step, 0], v[step, 0], radiation)
        if step == len(time) - 1:
            break
        extrapolated = PchipInterpolator(
            radiation_times, radiation_values, extrapolate=True,
        )

        def derivative(stage_time, angle, speed):
            return speed, acceleration(
                stage_time, angle, speed, float(extrapolated(stage_time)),
            )

        k1 = derivative(at_time, q[step, 0], v[step, 0])
        k2 = derivative(at_time + dt / 2,
                        q[step, 0] + dt * k1[0] / 2,
                        v[step, 0] + dt * k1[1] / 2)
        k3 = derivative(at_time + dt / 2,
                        q[step, 0] + dt * k2[0] / 2,
                        v[step, 0] + dt * k2[1] / 2)
        k4 = derivative(at_time + dt,
                        q[step, 0] + dt * k3[0],
                        v[step, 0] + dt * k3[1])
        q[step + 1, 0] = q[step, 0] + dt * (
            k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]
        ) / 6
        v[step + 1, 0] = v[step, 0] + dt * (
            k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]
        ) / 6
        held.last_heading = current_heading

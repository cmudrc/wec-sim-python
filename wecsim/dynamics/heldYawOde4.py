"""Published WEC-Sim ode4 path for a single sampled passive-yaw coordinate.

The source samples convolution at major steps, extrapolates its last three
loads to RK4 minor stages. Held excitation delays its selected heading by one
major step; a nearest-heading bank selects at each RK4 stage. This opt-in path
preserves those numerical details without changing the implicit default.
"""

import numpy as np
from scipy.interpolate import PchipInterpolator

from wecsim.hydrodynamics.passiveYaw import (
    HeldPassiveYawExcitation, NearestSampledHeadingExcitation,
)


def integrate_held_yaw_ode4(system, time, q, v, a, dt):
    """Fill one-coordinate motion arrays using the source's sampled step."""
    moving = [index for index, body in enumerate(system.bodies)
              if isinstance(body.state_excitation, (
                  HeldPassiveYawExcitation, NearestSampledHeadingExcitation,
              ))]
    if (system.coordinate_count != 1 or len(moving) != 1
            or system.nonlinear_force is not None
            or system.radiation_discretization != "trapezoid"
            or any(body.state_inertia is not None
                   or body.radiation_kernel is None for body in system.bodies)):
        raise ValueError("wecsim_ode4 needs one pure sampled-yaw coordinate")
    moving_index = moving[0]
    selected = system.bodies[moving_index].state_excitation
    for index, body in enumerate(system.bodies):
        motion = body.motion(np.zeros(1), np.zeros(1))
        expected = np.zeros((6, 1))
        if index == moving_index:
            expected[5, 0] = 1
        if (not np.array_equal(motion.jacobian, expected)
                or not np.array_equal(motion.bias_acceleration, np.zeros(6))):
            raise ValueError("wecsim_ode4 needs pure yaw and stationary other bodies")
    if (len(selected.model.time) != 2 * len(time) - 1
            or (len(time) > 1 and not np.isclose(
                selected.model.time[1], dt / 2, rtol=0, atol=1e-10))):
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

        # The held-heading source has a Memory block, so only that variant
        # commits its current selection after the RK4 stages. The direction
        # bank has no Memory block and selects from each stage's current yaw.
        if isinstance(selected, HeldPassiveYawExcitation):
            current_heading = selected._heading(float(q[step, 0]))
            selected.force_history.append(selected._force(
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
        if isinstance(selected, HeldPassiveYawExcitation):
            selected.last_heading = current_heading

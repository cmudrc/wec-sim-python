"""Source-style adaptive regular-wave motion with delayed added mass.

The optional major-step schedule is a diagnostic: it supplies MATLAB's
accepted times, but never its body states or forces. A replay using that
schedule does not establish independent published-case parity.
"""

import numpy as np
from scipy.integrate import RK45


class _PrescribedRK45(RK45):
    def _estimate_error_norm(self, K, h, scale):
        # Source times are already fixed by the diagnostic input. Retain the
        # RK45 stages and bypass error-based rejection for these intervals.
        return 0.0


def integrate_source_regular_ode45(device, time, coordinate, speed,
                                   acceleration, dt, major_times=None):
    """Advance regular-wave motion and extrapolate delayed mass at major steps."""
    n = device.coordinate_count
    zero_delay = tuple(np.zeros(6) for _ in device.bodies)
    device.delayed_body_acceleration = zero_delay
    acceleration[0] = device.acceleration(time[0], coordinate[0], speed[0])
    accepted_time = [float(time[0])]

    def body_acceleration(q, v, qddot):
        result = []
        for body in device.bodies:
            motion = body.motion(q, v)
            result.append(motion.jacobian @ qddot + motion.bias_acceleration)
        return np.stack(result)

    accepted_body_acceleration = [body_acceleration(
        coordinate[0], speed[0], acceleration[0],
    )]

    def derivative(at_time, state):
        if at_time <= accepted_time[-1] + 1e-14:
            delayed = (zero_delay if len(accepted_time) == 1
                       else device.delayed_body_acceleration)
        elif len(accepted_time) == 1:
            delayed = accepted_body_acceleration[-1]
        else:
            previous_interval = accepted_time[-1] - accepted_time[-2]
            fraction = ((at_time - accepted_time[-1] - device.added_mass_delay)
                        / previous_interval)
            delayed = (accepted_body_acceleration[-1] + fraction *
                       (accepted_body_acceleration[-1]
                        - accepted_body_acceleration[-2]))
        device.delayed_body_acceleration = tuple(delayed)
        q, v = state[:n], state[n:]
        try:
            return np.r_[v, device.acceleration(at_time, q, v)]
        except ValueError as exc:
            raise RuntimeError(f"source-style regular dynamics failed at {at_time:.9g} s") from exc

    if len(time) == 1:
        return
    if major_times is None:
        targets = time[1:]
    else:
        major_times = np.asarray(major_times, dtype=float)
        if (major_times.ndim != 1 or len(major_times) < 2
                or not np.isfinite(major_times).all()
                or not np.all(np.diff(major_times) > 0)
                or abs(major_times[0] - time[0]) > 1e-10
                or abs(major_times[-1] - time[-1]) > 1e-9):
            raise ValueError("major steps must increase from start to end time")
        targets = major_times[1:]

    def save_sample(sample, target, state, slope):
        if sample < len(time) and abs(target - time[sample]) < 1e-9:
            coordinate[sample], speed[sample], acceleration[sample] = (
                state[:n], state[n:], slope[n:]
            )
            return sample + 1
        return sample

    solver = _PrescribedRK45 if major_times is not None else RK45
    stepper = solver(
        derivative, float(time[0]), np.r_[coordinate[0], speed[0]],
        float(time[-1]), first_step=min(dt, float(targets[0] - time[0])),
        max_step=dt, rtol=1e-3, atol=1e-6,
    )
    sample = 1
    for target in targets:
        target = float(target)
        while stepper.t < target - 1e-12:
            stepper.max_step = target - stepper.t
            if major_times is not None:
                stepper.h_abs = target - stepper.t
            stepper.step()
            if stepper.status == "failed":
                raise RuntimeError("source-style adaptive regular solver failed")
            if major_times is not None and stepper.t < target - 1e-12:
                raise RuntimeError("a prescribed major step was rejected")
            accepted_time.append(float(stepper.t))
            accepted_body_acceleration.append(body_acceleration(
                stepper.y[:n], stepper.y[n:], stepper.f[n:],
            ))
        if abs(stepper.t - target) > 1e-9:
            raise RuntimeError("adaptive solver missed a requested step")
        sample = save_sample(sample, target, stepper.y, stepper.f)
    if sample != len(time):
        raise ValueError("major-step schedule must contain every output time")

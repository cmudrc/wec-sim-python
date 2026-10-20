"""Source-style adaptive regular-wave motion with delayed added mass.

The optional major-step schedule is a diagnostic: it supplies MATLAB's
accepted times, but never its body states or forces. A replay using that
schedule does not establish independent published-case parity.
"""

import numpy as np
from scipy.integrate import RK45
from scipy.integrate._ivp.rk import rk_step


def _rm3_source_state(y):
    """Express the reduced RM3 coordinates in Simscape's eight logged states."""
    q0, q1, q2, angle, v0, v1, v2, spin = y
    sine, cosine = np.sin(angle), np.cos(angle)
    return np.array([
        q0 + q2 * sine, q2 * cosine, -angle,
        v0 + v2 * sine + q2 * cosine * spin,
        v2 * cosine - q2 * sine * spin, -spin,
        q1 - q2, v1 - v2,
    ])


class _SourceRK45(RK45):
    """Use source-state error control for independently chosen RM3 steps."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # With RelTol=1e-3, Simulink's auto AbsTol starts at 1e-6 and
        # increases with the largest magnitude each source state has reached.
        self.max_source_state = np.maximum(abs(_rm3_source_state(self.y)), 1e-3)

    def step(self):
        # MATLAB ode45 stretches a step within 10% of its target. Here the
        # target is the next published output time, which is a major step.
        if 1.1 * self.h_abs >= self.max_step:
            self.h_abs = self.max_step
        result = super().step()
        self.max_source_state = np.maximum(
            self.max_source_state, abs(_rm3_source_state(self.y)))
        return result

    def _step_impl(self):
        # SciPy supplies the Dormand-Prince stages; MATLAB ode45 uses a 0.8
        # safety factor and caps accepted growth at five, rather than SciPy's
        # 0.9 and ten. Keep these choices local to the source-style solver.
        t, y = self.t, self.y
        min_step = 10 * abs(np.nextafter(t, self.direction * np.inf) - t)
        h_abs = min(self.max_step, max(min_step, self.h_abs))
        rejected = False
        while True:
            if h_abs < min_step:
                return False, self.TOO_SMALL_STEP
            t_new = t + h_abs * self.direction
            if self.direction * (t_new - self.t_bound) > 0:
                t_new = self.t_bound
            h = t_new - t
            h_abs = abs(h)
            y_new, f_new = rk_step(
                self.fun, t, y, self.f, h, self.A, self.B, self.C, self.K,
            )
            error = self._estimate_error_norm(self.K, h, None)
            if error < 1:
                factor = (5.0 if error == 0 else
                          min(5.0, 0.8 * error ** self.error_exponent))
                if rejected:
                    factor = min(1.0, factor)
                h_abs *= factor
                break
            h_abs *= max(0.1, 0.8 * error ** self.error_exponent)
            rejected = True
        self.h_previous = h
        self.y_old = y
        self.t, self.y = t_new, y_new
        self.h_abs, self.f = h_abs, f_new
        return True, None

    def _estimate_error_norm(self, K, h, scale):
        y5 = self.y + h * (K[:-1].T @ self.B)
        y4 = y5 - self._estimate_error(K, h)
        x5, x4 = _rm3_source_state(y5), _rm3_source_state(y4)
        tolerance = self.rtol * np.maximum.reduce((
            abs(x5), abs(x4), self.max_source_state,
        ))
        return float(np.max(abs(x5 - x4) / tolerance))


class _PrescribedRK45(RK45):
    def _estimate_error_norm(self, K, h, scale):
        # Source times are already fixed by the diagnostic input. Retain the
        # RK45 stages and bypass error-based rejection for these intervals.
        return 0.0


def integrate_source_regular_ode45(device, time, coordinate, speed,
                                   acceleration, dt, major_times=None):
    """Advance regular-wave motion and extrapolate delayed mass at major steps."""
    n = device.coordinate_count
    if n != 4 or len(device.bodies) != 2:
        raise ValueError("source regular ode45 needs the two-body RM3 joint")
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

    solver = _PrescribedRK45 if major_times is not None else _SourceRK45
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

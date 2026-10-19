"""Actuator used by the published WaveStar nonlinear predictive controller."""

import numpy as np
from scipy.optimize import lsq_linear
from scipy.signal import cont2discrete


def _discrete_plant(dt):
    """Return the published controller's four-state model and moment input."""
    inertia = 1.04 + .4805
    continuous = np.array([
        [0., 1., 0., 0.],
        [-92.33 / inertia, -(-.1586 + 1.8) / inertia,
         -4.739 / inertia, -.5 / inertia],
        [0., 8., -13.59, -13.35],
        [0., 0., 8., 0.],
    ])
    input_matrix = np.array([[0.], [1 / inertia], [0.], [0.]])
    transition, input_discrete, _, _, _ = cont2discrete(
        (continuous, input_matrix, np.eye(4), np.zeros((4, 1))),
        dt, method="zoh",
    )
    return transition, input_discrete[:, 0]


class WaveStarNmpcActuator:
    """Map a sampled torque request to axial PTO force.

    The controller's kinematic constants are retained separately from the
    published body's B-to-C geometry. The Simulink model uses these values in
    its torque-to-force block, even though its physical rod is slightly
    shorter at the neutral pose.
    """

    def __init__(self):
        self._command = [0.0, 0.0]
        self._torque = [0.0, 0.0]

    @staticmethod
    def moment_arm(stroke):
        """Return the controller's idealized metres-per-radian arm."""
        length = np.asarray(stroke, dtype=float) + .381408
        cosine = (length**2 - .412**2 - .2**2) / (-2 * .412 * .2)
        if (not np.isfinite(length).all() or np.any(length <= 0)
                or np.any(np.abs(cosine) >= 1)):
            raise ValueError("WaveStar NMPC stroke is outside the controller linkage")
        return .412 * .2 * np.sqrt(1 - cosine**2) / length

    def step(self, command_torque: float, stroke: float) -> float:
        """Advance one published 0.05 s actuator sample and return force N."""
        if not np.isfinite(command_torque):
            raise ValueError("WaveStar NMPC torque command must be finite")
        arm = self.moment_arm(stroke)
        command = float(np.clip(command_torque, -12, 12))
        torque = (
            .96294775 * command
            - 1.92342278 * self._command[0]
            + .96053421 * self._command[1]
            + 1.99777700 * self._torque[0]
            - .99783206 * self._torque[1]
        )
        self._command[:] = [command, self._command[0]]
        self._torque[:] = [torque, self._torque[0]]
        return float(np.clip(torque, -12, 12) / arm)

    def force_at(self, stroke: float) -> float:
        """Project the held motor torque through the current rod geometry."""
        return float(np.clip(self._torque[0], -12, 12)
                     / self.moment_arm(stroke))


class WaveStarNmpcObserver:
    """Replay the published sampled position/velocity filter and Kalman observer.

    ``step`` takes the measured PTO stroke and the preceding torque command.
    It returns pitch, pitch speed, two radiation states, and estimated wave
    excitation moment. The observer uses the controller's idealized linkage,
    which is slightly different from the physical B-to-C rod geometry.
    """

    def __init__(self, dt: float = .05):
        if not np.isfinite(dt) or dt <= 0:
            raise ValueError("WaveStar NMPC observer step must be positive")
        self.dt = float(dt)
        transition, input_discrete = _discrete_plant(self.dt)
        self._transition = np.eye(5)
        self._transition[:4, :4] = transition
        self._transition[:4, 4] = input_discrete
        self._input = np.r_[input_discrete, 0.]
        self._measurement = np.eye(5)[:2]
        self._process_covariance = np.diag([5e-6] * 4 + [.35])
        self._measurement_covariance = np.diag([6e-6] * 2)
        self._covariance = self._process_covariance.copy()
        self._state = np.zeros(5)
        self._previous_position = 0.
        self._previous_highpass = 0.
        self._previous_velocity = 0.

    @staticmethod
    def _position_from_stroke(stroke: float) -> float:
        length = float(stroke) + .381408
        cosine = (length**2 - .412**2 - .2**2) / (-2 * .412 * .2)
        if (not np.isfinite(length) or length <= 0
                or not np.isfinite(cosine) or abs(cosine) >= 1):
            raise ValueError("WaveStar NMPC stroke is outside the controller linkage")
        return float(np.arccos(cosine) - 1.170165)

    def step(self, stroke: float, previous_command: float) -> np.ndarray:
        """Advance one sample; ``previous_command`` is torque in N m."""
        if not np.isfinite(previous_command):
            raise ValueError("WaveStar NMPC previous command must be finite")
        position = self._position_from_stroke(stroke)
        frequency = 2 * np.pi * 30
        denominator = 2 + frequency * self.dt
        delayed_denominator = frequency * self.dt - 2
        highpass = (
            2 * frequency * (position - self._previous_position)
            - delayed_denominator * self._previous_highpass
        ) / denominator
        velocity = (
            frequency * self.dt * (highpass + self._previous_highpass)
            - delayed_denominator * self._previous_velocity
        ) / denominator
        self._previous_position = position
        self._previous_highpass = highpass
        self._previous_velocity = velocity

        predicted = (self._transition @ self._state
                     + self._input * previous_command)
        covariance = (self._transition @ self._covariance
                      @ self._transition.T + self._process_covariance)
        innovation_covariance = (
            self._measurement @ covariance.T @ self._measurement.T
            + self._measurement_covariance
        )
        gain = np.linalg.solve(
            innovation_covariance, self._measurement @ covariance.T,
        ).T
        measurement = np.array([position, velocity])
        self._state = (predicted + gain @ (
            measurement - self._measurement @ predicted
        ))
        self._covariance = covariance - gain @ self._measurement @ covariance
        # The pinned observer applies this calibration to its persistent
        # excitation state after every update, including the first sample.
        self._state[-1] -= .6979
        return self._state.copy()


class WaveStarNmpcPredictor:
    """The published forward-backward AR forecast of excitation moment.

    The source keeps ``training_set * order`` past samples, refits the AR
    coefficients at each integer second, and forecasts ``horizon`` samples.
    Defaults reproduce Sea State 6 of the pinned WaveStar NMPC application.
    """

    def __init__(
        self, *, dt: float = .05, order: int = 18,
        training_set: int = 10, horizon: int = 40,
        start_time: float = 10.,
    ):
        if (not np.isfinite([dt, start_time]).all() or dt <= 0
                or start_time < 0 or not all(
                    isinstance(value, int) and value > 0
                    for value in (order, training_set, horizon)
                )):
            raise ValueError("WaveStar NMPC predictor settings are invalid")
        window = order * training_set
        samples_per_second = round(1 / dt)
        start_sample = round(start_time / dt)
        if (window <= order or samples_per_second < 1
                or not np.isclose(samples_per_second * dt, 1,
                                  rtol=0, atol=1e-12)
                or not np.isclose(start_sample * dt, start_time,
                                  rtol=0, atol=1e-12)):
            raise ValueError("WaveStar NMPC predictor needs aligned samples")
        self.dt = float(dt)
        self.order = order
        self.horizon = horizon
        self._samples_per_second = samples_per_second
        self._start_sample = start_sample
        self._past = np.zeros(window)
        self._coefficients = None
        self._sample = -1

    def step(self, estimated_excitation_moment: float) -> np.ndarray:
        """Advance one sample and return the future moment, in N m."""
        if not np.isfinite(estimated_excitation_moment):
            raise ValueError("WaveStar NMPC estimated moment must be finite")
        self._sample += 1
        self._past[:-1] = self._past[1:]
        self._past[-1] = estimated_excitation_moment
        forecast = np.zeros(self.horizon)
        if self._sample < self._start_sample:
            return forecast

        if (self._coefficients is None
                or self._sample % self._samples_per_second == 0):
            # MATLAB ar(y, order) defaults to the modified-covariance,
            # forward-backward fit without windowing. Both regressions use
            # the same coefficient vector, with the second series reversed.
            p = self.order
            n = len(self._past)
            forward = np.array([
                self._past[t - p:t][::-1] for t in range(p, n)
            ])
            backward = np.array([
                self._past[t - p + 1:t + 1] for t in range(p, n)
            ])
            targets = np.r_[self._past[p:], self._past[:n - p]]
            self._coefficients = np.linalg.lstsq(
                np.vstack((forward, backward)), targets, rcond=None,
            )[0]

        history = list(self._past[-self.order:])
        for index in range(self.horizon):
            prediction = self._coefficients @ np.asarray(
                history[-self.order:][::-1]
            )
            forecast[index] = prediction
            history.append(prediction)
        return forecast


class WaveStarNmpcController:
    """Replay the published Sea State 6 resistive startup and RTI NMPC.

    ``step`` takes the five observer states and an excitation forecast for
    the current sample. It returns the requested PTO torque in N m. The
    controller shifts its prior horizon solution before each QP, matching
    the source's persistent nominal input. It does not integrate the WEC.
    """

    def __init__(
        self, *, dt: float = .05, horizon: int = 40,
        gain: float = 19.4, torque_limit: float = 12.,
        input_penalty: float = .9, generator_efficiency: float = .7,
        motoring_efficiency: float = 1 / .7, smoothness: float = 1000.,
        start_time: float = 15.,
    ):
        parameters = [dt, gain, torque_limit, input_penalty,
                      generator_efficiency, motoring_efficiency,
                      smoothness, start_time]
        if (not np.isfinite(parameters).all() or dt <= 0
                or not isinstance(horizon, int) or horizon < 1
                or gain < 0 or torque_limit <= 0 or input_penalty <= 0
                or generator_efficiency <= 0 or motoring_efficiency <= 0
                or smoothness <= 0 or start_time < 0):
            raise ValueError("WaveStar NMPC controller settings are invalid")
        self.dt = float(dt)
        self.horizon = horizon
        self.gain = gain
        self.torque_limit = torque_limit
        self.start_time = start_time
        self._alpha = (motoring_efficiency + generator_efficiency) / 2
        self._beta = (motoring_efficiency - generator_efficiency) / 2
        self._smoothness = smoothness
        self._weight = np.kron(np.eye(horizon), [[0., 1.], [1., 0.]])
        self._penalty = input_penalty * np.eye(horizon)
        transition, input_discrete = _discrete_plant(dt)
        self._transition = np.zeros((5, 5))
        self._transition[:4, :4] = transition
        self._input = np.r_[input_discrete, 1.]
        self._wave_input = np.r_[input_discrete, 0.]
        self._nominal = np.zeros(horizon)

    def step(self, estimated_state, excitation_forecast, time: float) -> float:
        """Return one sampled command; ``time`` is seconds from simulation start."""
        state = np.asarray(estimated_state, dtype=float)
        forecast = np.asarray(excitation_forecast, dtype=float)
        if (state.shape != (5,) or forecast.shape != (self.horizon,)
                or not np.isfinite(state).all()
                or not np.isfinite(forecast).all()
                or not np.isfinite(time) or time < 0):
            raise ValueError("WaveStar NMPC controller needs finite sampled inputs")
        if 10 <= time < self.start_time:
            return float(np.clip(-self.gain * state[1],
                                 -self.torque_limit, self.torque_limit))
        if time < self.start_time:
            return 0.

        nominal = np.r_[self._nominal[1:], self._nominal[-1]]
        state = state.copy()
        sensitivity = np.zeros((5, self.horizon))
        output_sensitivity = np.zeros((2 * self.horizon, self.horizon))
        output = np.zeros(2 * self.horizon)
        for index in range(self.horizon):
            state = (self._transition @ state + self._input * nominal[index]
                     + self._wave_input * forecast[index])
            speed, command = state[1], state[4]
            product = np.tanh(self._smoothness * command * speed)
            efficiency = self._alpha + self._beta * product
            gradient = np.zeros((2, 5))
            gradient[0, 1] = 1.
            gradient[1, 1] = (self._smoothness * self._beta * command**2
                              * (1 - product**2))
            gradient[1, 4] = (efficiency + self._smoothness * self._beta
                              * command * speed * (1 - product**2))
            sensitivity = self._transition @ sensitivity
            sensitivity[:, index] += self._input
            row = slice(2 * index, 2 * index + 2)
            output_sensitivity[row] = gradient @ sensitivity
            output[row] = [speed, command * efficiency]

        hessian = (self._penalty
                   + output_sensitivity.T @ self._weight @ output_sensitivity)
        hessian = (hessian + hessian.T) / 2
        linear = output_sensitivity.T @ self._weight @ output
        # H is positive definite because the input penalty is positive.
        # Factor the quadratic into a bounded least-squares problem so the
        # active-set solve reaches its constrained optimum directly.
        factor = np.linalg.cholesky(hessian).T
        target = -np.linalg.solve(factor.T, linear)
        lower = -self.torque_limit - nominal
        upper = self.torque_limit - nominal
        solution = lsq_linear(
            factor, target, bounds=(lower, upper), method="bvls",
            tol=1e-12, max_iter=100,
        )
        if not solution.success:
            raise RuntimeError(f"WaveStar NMPC QP failed: {solution.message}")
        self._nominal = nominal + solution.x
        return float(np.clip(self._nominal[0],
                             -self.torque_limit, self.torque_limit))


class WaveStarNmpcPTO:
    """Sampled WaveStar NMPC PTO for ``run_wavestar_published``.

    The plant advances at ``plant_dt`` while the observer, predictor,
    controller, and actuator update every ``control_dt``. The resulting
    motor torque is held between controller samples while the axial force
    follows the moving rod moment arm. Defaults reproduce the published
    Sea State 6 control settings with a fine plant step.
    """

    def __init__(
        self, *, plant_dt: float = .001, control_dt: float = .05,
        observer: WaveStarNmpcObserver | None = None,
        predictor: WaveStarNmpcPredictor | None = None,
        controller: WaveStarNmpcController | None = None,
        actuator: WaveStarNmpcActuator | None = None,
    ):
        if (not np.isfinite([plant_dt, control_dt]).all()
                or plant_dt <= 0 or control_dt <= 0):
            raise ValueError("WaveStar NMPC sample steps must be positive")
        ratio = round(control_dt / plant_dt)
        if (ratio < 1 or not np.isclose(ratio * plant_dt, control_dt,
                                       rtol=0, atol=1e-12)):
            raise ValueError("WaveStar NMPC control step must align with plant")
        for component in (observer, predictor, controller):
            if (component is not None and not np.isclose(
                    component.dt, control_dt, rtol=0, atol=1e-12)):
                raise ValueError("WaveStar NMPC component step must match control")
        self.plant_dt = float(plant_dt)
        self.control_dt = float(control_dt)
        self.observer = observer or WaveStarNmpcObserver(dt=control_dt)
        self.predictor = predictor or WaveStarNmpcPredictor(dt=control_dt)
        self.controller = controller or WaveStarNmpcController(dt=control_dt)
        self.actuator = actuator or WaveStarNmpcActuator()
        if self.predictor.horizon != self.controller.horizon:
            raise ValueError("WaveStar NMPC forecast and control horizons differ")
        self._ratio = ratio
        self._step = 0
        self._force = 0.
        self._previous_command = 0.
        self._command_time = []
        self._command_torque = []

    @property
    def command_time(self) -> np.ndarray:
        """Times at which the PTO command was updated, in seconds."""
        return np.asarray(self._command_time)

    @property
    def command_torque(self) -> np.ndarray:
        """Requested PTO torque at each controller sample, in N m."""
        return np.asarray(self._command_torque)

    def step(self, stroke: float, *, noise: float = 0.,
             dropout: bool = False) -> float:
        """Advance one plant sample and return the held axial PTO force."""
        if noise != 0 or dropout:
            raise ValueError("WaveStar NMPC sensor faults are unsupported")
        if self._step % self._ratio == 0:
            instant = self._step * self.plant_dt
            state = self.observer.step(stroke, self._previous_command)
            forecast = self.predictor.step(state[-1])
            command = self.controller.step(state, forecast, instant)
            self._force = self.actuator.step(command, stroke)
            self._previous_command = command
            self._command_time.append(instant)
            self._command_torque.append(command)
        else:
            self._force = self.actuator.force_at(stroke)
        self._step += 1
        return self._force

    def force_at_stroke(self, stroke: float) -> float:
        """Read the current torque at an intermediate plant stroke."""
        return self.actuator.force_at(stroke)

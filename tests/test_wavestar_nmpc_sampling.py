"""Portable guards for the published WaveStar controller's sample timing."""

import numpy as np

from wecsim import WaveStarNmpcController, WaveStarNmpcPTO


def test_controller_startup_switches_from_idle_to_resistive_control():
    controller = WaveStarNmpcController(horizon=2)
    forecast = np.zeros(2)
    state = np.array([0., 1., 0., 0., 0.])

    assert controller.step(state, forecast, 9.95) == 0.
    assert controller.step(state, forecast, 10.) == -12.
    state[1] = -.25
    assert controller.step(state, forecast, 14.95) == 4.85


class _Observer:
    dt = .05

    def __init__(self):
        self.previous_commands = []

    def step(self, stroke, previous_command):
        self.previous_commands.append(previous_command)
        return np.array([0., 0., 0., 0., stroke])


class _Predictor:
    dt = .05
    horizon = 2

    def step(self, estimated_moment):
        return np.full(self.horizon, estimated_moment)


class _Controller:
    dt = .05
    horizon = 2

    def __init__(self):
        self.times = []

    def step(self, state, forecast, time):
        self.times.append(time)
        np.testing.assert_allclose(forecast, state[-1])
        return 3. + 2 * (len(self.times) - 1)


def test_pto_hands_previous_command_to_observer_at_control_samples():
    observer = _Observer()
    controller = _Controller()
    pto = WaveStarNmpcPTO(
        plant_dt=.01, control_dt=.05, observer=observer,
        predictor=_Predictor(), controller=controller,
    )

    force = [pto.step(.01 * sample) for sample in range(11)]

    np.testing.assert_allclose(pto.command_time, [0., .05, .1], atol=1e-15)
    np.testing.assert_array_equal(pto.command_torque, [3., 5., 7.])
    np.testing.assert_array_equal(observer.previous_commands, [0., 3., 5.])
    np.testing.assert_allclose(controller.times, pto.command_time)
    assert np.isfinite(force).all()

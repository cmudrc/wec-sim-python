"""The WaveStar NMPC actuator holds torque, not axial force."""

import numpy as np

from wecsim import WaveStarNmpcActuator, WaveStarNmpcPTO


class _Observer:
    dt = .05

    def step(self, stroke, previous_command):
        return np.zeros(5)


class _Predictor:
    dt = .05
    horizon = 40

    def step(self, estimated_moment):
        return np.zeros(self.horizon)


class _Controller:
    dt = .05
    horizon = 40

    def step(self, state, forecast, time):
        return 6.


def test_held_motor_torque_tracks_intersample_rod_stroke():
    pto = WaveStarNmpcPTO(
        plant_dt=.001, control_dt=.05,
        observer=_Observer(), predictor=_Predictor(),
        controller=_Controller(),
    )
    first_stroke, second_stroke = 0., .02
    first_force = pto.step(first_stroke)
    second_force = pto.step(second_stroke)
    first_arm = WaveStarNmpcActuator.moment_arm(first_stroke)
    second_arm = WaveStarNmpcActuator.moment_arm(second_stroke)
    assert pto.command_time.tolist() == [0.]
    assert pto.command_torque.tolist() == [6.]
    assert abs(second_force - first_force) > .01
    np.testing.assert_allclose(
        second_force, first_force * first_arm / second_arm,
        rtol=0, atol=1e-12,
    )

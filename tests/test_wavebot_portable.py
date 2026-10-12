"""A short prescribed-gain WaveBot run using portable hydrodynamics."""

from pathlib import Path

import numpy as np
from scipy.io import savemat

from wecsim.control.wavebotControl import (
    WaveBotPIController, run_wavebot_prescribed_gains,
)
from wecsim.control.wavebotImpedance import run_wavebot_impedance
from wecsim.environment.irregularWave import pm_equal_energy_components


SPHERE = Path(__file__).resolve().parents[1] / "examples/data/sphere.h5"


def test_prescribed_wavebot_gain_history_replays_actuator_force():
    sea = pm_equal_energy_components(
        SPHERE, significant_height=1, peak_period=8,
        directions=[0], spreading=[1], count=8, seed=7,
    )
    gains = np.tile([-1000., 0., -1000., 500., -200., 0.], (21, 1))
    result = run_wavebot_prescribed_gains(
        SPHERE, sea, gains, dt=.01, end_time=.2, ramp_time=0,
    )
    assert result.body_position.shape == result.body_velocity.shape == (21, 6)
    assert result.actuator_command.shape == (21, 3)
    assert np.isfinite(result.body_position).all()
    assert np.isfinite(result.actuator_command).all()
    np.testing.assert_array_equal(result.gains, gains)
    controller = WaveBotPIController()
    for index, speed in enumerate(result.body_velocity[:, [0, 2, 4]]):
        controller.set_gains(gains[index])
        np.testing.assert_allclose(result.actuator_command[index],
                                   controller.force(speed), rtol=0, atol=1e-9)
        controller.advance(speed, .01)


def test_impedance_runner_scales_multisine_and_keeps_damping_dissipative(tmp_path):
    source = tmp_path / "multisine.mat"
    savemat(source, {"multisine": {
        "time": np.array([0., .1, .2]),
        "signals": {"values": np.array([
            [0., 0., 0.], [1., 2., 3.], [2., 4., 6.],
        ])},
    }})
    result = run_wavebot_impedance(
        SPHERE, source, dt=.01, end_time=.2,
    )
    assert result.body_position.shape == (21, 6)
    np.testing.assert_allclose(result.applied_force[10], [50, 200, -12])
    np.testing.assert_allclose(result.applied_force[20], [100, 400, -24])
    assert np.isfinite(result.body_position).all()
    damping_power = np.sum(result.linear_damping_force * result.body_velocity,
                           axis=1)
    assert np.max(damping_power) <= 1e-8
    np.testing.assert_allclose(result.mooring_force[:, 0],
                               -24000 * result.body_position[:, 0]
                               - 5 * result.body_velocity[:, 0],
                               rtol=0, atol=1e-8)

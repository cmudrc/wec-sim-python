"""Portable controller checks that do not require a MATLAB trajectory."""

import numpy as np
import pytest

from wecsim.control.mostController import MostBaselineController
from wecsim.control.wavebotControl import WaveBotPIController, WaveBotPowerAutotuner


def test_wavebot_pi_uses_measured_speed_then_advances_integral():
    control = WaveBotPIController((2, 3, -4, 5, 6, -7))
    speed = np.array([1., -2., 3.])
    np.testing.assert_array_equal(control.force(speed), [2, 8, 18])
    control.advance(speed, .25)
    np.testing.assert_array_equal(control.force(speed), [2.75, 5.5, 12.75])
    with pytest.raises(ValueError, match="positive dt"):
        control.advance(speed, 0)
    with pytest.raises(ValueError, match="six finite"):
        control.set_gains([1, 2])


def test_wavebot_excitation_estimate_subtracts_delayed_command():
    frequency = np.array([.2, 1.])
    impedance = np.zeros((3, 3, 2), dtype=complex)
    impedance[0, 0] = impedance[1, 1] = impedance[2, 2] = 4 + 2j
    tuner = WaveBotPowerAutotuner(impedance, frequency)
    speed = np.zeros((1024, 3))
    command = np.zeros_like(speed)
    speed[-1, 0] = 2
    command[-1, 0] = 3
    # The final Hamming sample weights the pulse; the FFT of a pulse has
    # constant magnitude and a frequency-dependent phase.
    pulse = np.fft.fft(np.hamming(1024) * np.r_[np.zeros(1023), 1])[52:257]
    expected = 2 * (2 * (4 + 2j) - 3) * pulse / 1024
    np.testing.assert_allclose(tuner.excitation_spectrum(speed, command)[:, 0],
                               expected, rtol=0, atol=1e-14)
    np.testing.assert_array_equal(tuner.excitation_spectrum(speed, command)[:, 1:], 0)
    with pytest.raises(ValueError, match="matching finite"):
        tuner.excitation_spectrum(speed, command[:-1])


def test_most_controller_obeys_torque_and_pitch_rate_limits():
    control = MostBaselineController.iea15mw()
    time = np.arange(21) * .05
    rotor_speed = np.r_[control.initial_omega, np.full(20, 1.2)]
    torque, pitch = control.simulate(time, rotor_speed)
    assert np.isfinite(torque).all() and np.isfinite(pitch).all()
    assert np.max(np.abs(np.diff(torque))) <= control.torque_max_rate * .05 + 1e-7
    assert np.max(np.abs(np.diff(pitch))) <= control.pitch_max_rate * .05 + 1e-12
    assert np.all((pitch >= 0) & (pitch <= np.pi / 2))
    assert torque[-1] > torque[0]
    assert pitch[-1] > pitch[0]
    with pytest.raises(ValueError, match="aligned"):
        control.simulate([0, .1], [1])

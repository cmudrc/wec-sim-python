"""Keep spatial irregular-wave gauges paired to a saved 200 s history."""

from pathlib import Path

import numpy as np

from wecsim.environment.waveClass import WaveClass


FIXTURE = (Path(__file__).resolve().parent
           / "test_objects/test_waveclass/testData/irregular_1_test")


def test_equal_energy_irregular_wave_matches_three_spatial_gauges():
    wave = WaveClass("irregular")
    wave.freqDisc = "EqualEnergy"
    wave.numFreq = 500
    wave.T = 8
    wave.H = 2.5
    wave.spectrumType = "BS"
    wave.waveDir = [0]
    wave.waveSpread = [1]
    # Replay the saved realization so this tests wave and gauge synthesis,
    # independently of the random-number generator used to obtain phases.
    wave.phaseData = np.loadtxt(FIXTURE / "phase.txt")[:, None]
    wave.wavegauge1loc = [5, 5]
    wave.wavegauge2loc = [10, 0]
    wave.wavegauge3loc = [0, -10]
    wave.waveSetup(
        [5.19999512307279, 0.0199999977946844], "infinite",
        100, .1, 2000, 9.81, 1000, 200,
    )

    for name in ("waveAmpTime", "waveAmpTime1", "waveAmpTime2",
                 "waveAmpTime3"):
        expected = np.loadtxt(FIXTURE / f"{name}.txt").T
        actual = np.asarray(getattr(wave, name))
        assert actual.shape == expected.shape == (2, 2001)
        np.testing.assert_allclose(actual, expected, rtol=0, atol=5e-12)
    assert np.max(np.abs(np.asarray(wave.waveAmpTime1[1])
                         - np.asarray(wave.waveAmpTime2[1]))) > .01

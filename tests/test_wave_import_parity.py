"""Pair production wave imports with archived MATLAB-generated histories."""

from pathlib import Path

import numpy as np
import pytest

from wecsim.environment.waveClass import WaveClass


FIXTURES = (Path(__file__).resolve().parent / "test_objects"
            / "test_waveclass" / "testData")
BEM_FREQUENCIES = [5.19999512307279, 0.0199999977946844]


@pytest.mark.parametrize("extension", ("txt", "mat"))
def test_imported_elevation_replays_matlab_interpolation_and_ramp(extension):
    fixture = FIXTURES / "etaImport_1_test"
    wave = WaveClass("etaImport")
    wave.etaDataFile = str(fixture / f"etaData.{extension}")
    wave.waveSetup(BEM_FREQUENCIES, "infinite", 100, .1, 4000,
                   9.81, 1000, 400)

    expected = np.loadtxt(fixture / "waveAmpTime.txt").T
    expected_gauge = np.loadtxt(fixture / "waveAmpTime1.txt").T
    np.testing.assert_allclose(wave.waveAmpTime, expected, rtol=0, atol=2e-12)
    for history in (wave.waveAmpTime1, wave.waveAmpTime2,
                    wave.waveAmpTime3):
        np.testing.assert_allclose(history, expected_gauge,
                                   rtol=0, atol=1e-12)
    assert wave.waveAmpTime[1][0] == 0
    assert np.max(np.abs(wave.waveAmpTime[1][1000:])) > .1


@pytest.mark.parametrize("extension", ("txt", "mat"))
def test_imported_spectrum_replays_matlab_bins_phase_and_sea(extension):
    fixture = FIXTURES / "spectrumImport_1_test"
    wave = WaveClass("spectrumImport")
    wave.spectrumDataFile = str(fixture / f"spectrumData.{extension}")
    wave.waveSetup(BEM_FREQUENCIES, "infinite", 100, .1, 2000,
                   9.81, 1000, 200)

    assert wave.numFreq == 1001
    for field in ("w", "k", "A", "S", "dw", "phase"):
        expected = np.loadtxt(fixture / f"{field}.txt")
        np.testing.assert_allclose(np.asarray(getattr(wave, field)).ravel(),
                                   expected, rtol=0, atol=1e-12)
    expected_sea = np.loadtxt(fixture / "waveAmpTime.txt").T
    np.testing.assert_allclose(wave.waveAmpTime, expected_sea,
                               rtol=0, atol=1e-12)

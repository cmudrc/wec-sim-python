"""Compare production Python output with a pinned, executed MATLAB waveClass."""

import os
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]

from wecsim.environment.waveClass import WaveClass  # noqa: E402
from wecsim.environment.irregularWave import (  # noqa: E402
    jonswap_equal_energy_components, pm_equal_energy_components,
    synthesize_irregular_response,
)

REFERENCE = os.environ.get("WEC_SIM_MATLAB_REFERENCE_DIR")
PORTABLE = ROOT / "tests/fixtures/matlab/wave_reference.zip"


def _source(name: str, *, ndmin: int = 0):
    """Read fresh MATLAB output when available, otherwise the pinned archive."""
    if REFERENCE:
        return np.loadtxt(Path(REFERENCE) / name, delimiter=",", ndmin=ndmin)
    with ZipFile(PORTABLE) as archive, archive.open(name) as values:
        return np.loadtxt(values, delimiter=",", ndmin=ndmin)


def test_regular_wave_against_executed_matlab():
    wave = WaveClass("regular")
    wave.T = 8.0
    wave.H = 2.5
    wave.waveDir = [0]
    wave.wavegauge1loc = [5, 5]
    wave.wavegauge2loc = [10, 0]
    wave.wavegauge3loc = [0, -10]
    wave.waveSetup(
        [5.19999512307279, 0.0199999977946844],
        "infinite", 100.0, 0.1, 2000, 9.81, 1000.0, 200.0,
    )
    origin = _source("regular_origin.csv")
    markers = _source("regular_markers.csv")
    np.testing.assert_allclose(np.asarray(wave.waveAmpTime).T, origin, rtol=0, atol=1e-10)
    for index, attribute in enumerate(("waveAmpTime1", "waveAmpTime2", "waveAmpTime3")):
        np.testing.assert_allclose(
            np.asarray(getattr(wave, attribute))[1], markers[:, index + 1],
            rtol=0, atol=1e-10,
        )


def test_finite_depth_against_executed_matlab():
    wave = WaveClass("regular")
    wave.T = 8.0
    wave.H = 2.5
    wave.waveDir = [0]
    wave.wavegauge1loc = [0, 0]
    wave.wavegauge2loc = [0, 0]
    wave.wavegauge3loc = [0, 0]
    wave.waveSetup([0.4, 2.0], 25.0, 0, 0.1, 2000, 9.81, 1000.0, 200.0)
    expected = _source("finite_depth.csv")
    np.testing.assert_allclose([np.asarray(wave.k).item(), np.asarray(wave.Pw).item()],
                               expected, rtol=1e-11)


@pytest.mark.parametrize("spectrum,height", [("PM", 2.5), ("JS", 4.0)])
@pytest.mark.parametrize(
    "discretization", ["Traditional", "EqualEnergy", "NarrowEqualEnergy"],
)
def test_current_irregular_spectrum_against_executed_matlab(
    spectrum, height, discretization,
):
    wave = WaveClass("irregular")
    wave.T = 8
    wave.H = height
    wave.spectrumType = spectrum
    wave.freqDisc = ("EqualEnergy" if discretization == "NarrowEqualEnergy"
                     else discretization)
    wave.numFreq = 64
    wave.phaseSeed = 1
    wave.phaseGenerator = "matlab"
    if spectrum == "PM":
        wave.waveDir = [0, 30, 90]
        wave.waveSpread = [0.1, 0.2, 0.7]
    wave.wavegauge1loc = [5, 5]
    wave.wavegauge2loc = [10, 0]
    wave.wavegauge3loc = [0, -10]
    suffix = {"Traditional": "", "EqualEnergy": "_equal",
              "NarrowEqualEnergy": "_narrow"}[discretization]
    label = spectrum.lower() + suffix
    if discretization == "NarrowEqualEnergy":
        wave.freqRange = [0.5, 1.5]
    wave.waveSetup([0.4, 2.0], "infinite", 1, 0.1, 20, 9.81, 1000, 2)
    source_phase = _source(f"{label}_phase.csv", ndmin=2)
    np.testing.assert_allclose(wave.phase.T, source_phase, rtol=0, atol=1e-13)

    if discretization != "Traditional":
        expected = _source(f"{label}_bins.csv")
        actual = np.column_stack((wave.w, wave.dw, wave.S))
    else:
        expected = _source(f"{label}_spectrum.csv")
        actual = np.column_stack((wave.w, wave.S))
    np.testing.assert_allclose(actual, expected,
                               rtol=2e-12, atol=1e-14)
    power = _source(f"{label}_power.csv")
    np.testing.assert_allclose(wave.Pw, power, rtol=2e-12, atol=1e-9)
    if spectrum == "JS":
        gamma = _source(f"{label}_gamma.csv")
        np.testing.assert_allclose(wave.gamma, gamma, rtol=2e-12)
    source_elevation = _source(f"{label}_elevation.csv")
    np.testing.assert_allclose(np.asarray(wave.waveAmpTime).T,
                               source_elevation, rtol=0, atol=2e-12)
    if discretization in ("Traditional", "NarrowEqualEnergy"):
        # The public WEC runner uses these component builders, while the
        # class above preserves Sungjun Won's original wave interface.
        hydro = ROOT / "tests/test_objects/test_bodyclass/testData/hydroData/oswec.h5"
        builder = (jonswap_equal_energy_components if spectrum == "JS"
                   else pm_equal_energy_components)
        components = builder(
            None if discretization == "Traditional" else hydro,
            significant_height=height, peak_period=8,
            directions=wave.waveDir, spreading=wave.waveSpread,
            count=64, phase=wave.phase.T,
            frequency_range=((0.4, 2.0) if discretization == "Traditional"
                             else (0.5, 1.5)),
            discretization=("traditional" if discretization == "Traditional"
                            else "equal_energy"),
        )
        actual_bins = np.column_stack((components.omega,
                                       components.spectral_amplitude / 2))
        expected_bins = expected
        if discretization != "Traditional":
            actual_bins = np.column_stack((components.omega, components.d_omega,
                                           components.spectral_amplitude / 2))
        np.testing.assert_allclose(
            actual_bins, expected_bins, rtol=2e-12, atol=1e-14,
        )
        if discretization == "Traditional":
            np.testing.assert_allclose(
                components.d_omega, np.full(64, 1.6 / 63), rtol=0, atol=1e-15,
            )
        incident = synthesize_irregular_response(
            hydro, components, dt=0.1, end_time=2, ramp_time=1,
        )
        np.testing.assert_allclose(
            incident.elevation, source_elevation[:, 1], rtol=0, atol=2e-12,
        )
    source_markers = _source(f"{label}_markers.csv")
    for index, attribute in enumerate(("waveAmpTime1", "waveAmpTime2",
                                       "waveAmpTime3")):
        np.testing.assert_allclose(
            np.asarray(getattr(wave, attribute))[1], source_markers[:, index + 1],
            rtol=0, atol=2e-12,
        )

"""Pair short RM3 wave-preprocessing branches with local pinned R2025b."""

from pathlib import Path

import numpy as np
import pytest

from wecsim import RegularCICWave, WEC
from wecsim.environment.waveClass import WaveClass


FIXTURE = Path(__file__).parent / "fixtures/matlab/core_cic_wave"


def _source(name: str) -> np.ndarray:
    return np.loadtxt(FIXTURE / name, delimiter=",", ndmin=2)


@pytest.mark.parametrize("kind", ("noWave", "noWaveCIC", "regularCIC"))
def test_core_wave_preprocessing_matches_pinned_matlab(kind):
    wave = WaveClass(kind)
    if kind != "noWaveCIC":
        wave.T = 8
    if kind == "regularCIC":
        wave.H = 2.5
        wave.waveDir = [30]
        wave.wavegauge1loc = [4, 0]
        wave.wavegauge2loc = [0, 5]
        wave.wavegauge3loc = [-3, 2]
    wave.waveSetup([0.4, 2], 25, 1, 0.1, 20, 9.81, 1000, 2)

    properties = np.array([wave.w, wave.T, wave.k, wave.A], dtype=float)
    np.testing.assert_allclose(properties, _source(f"{kind}_properties.csv").ravel(),
                               rtol=0, atol=1e-12)
    elevation = np.asarray(wave.waveAmpTime).T
    source_elevation = _source(f"{kind}_elevation.csv")
    assert elevation.shape == source_elevation.shape == (21, 2)
    np.testing.assert_allclose(elevation, source_elevation, rtol=0, atol=1e-12)

    if kind == "regularCIC":
        markers = np.column_stack((wave.waveAmpTime[0], *(getattr(
            wave, f"waveAmpTime{index}")[1] for index in (1, 2, 3))))
        np.testing.assert_allclose(markers, _source("regularCIC_markers.csv"),
                                   rtol=0, atol=1e-12)
        np.testing.assert_allclose(wave.Pw, _source("regularCIC_power.csv").item(),
                                   rtol=0, atol=1e-9)
        assert np.any(elevation[:, 1] != 0)
    else:
        assert np.all(elevation[:, 1] == 0)
        for index in (1, 2, 3):
            np.testing.assert_array_equal(
                getattr(wave, f"waveAmpTime{index}")[1], elevation[:, 1],
            )


def test_public_regular_cic_run_uses_source_incident_wave():
    """The installable builder must carry this short source sea to dynamics."""
    hydro = Path(__file__).resolve().parents[1] / "examples/data/oswec.h5"
    wec = WEC("regular CIC incident wave")
    flap = wec.body("flap", hydro)
    wec.coordinate("heave", flap.move("heave"))
    result = wec.run(
        RegularCICWave(2.5, 8, direction=30), dt=0.1,
        end_time=2, ramp_time=1, radiation_memory=0.2,
    )
    source = _source("regularCIC_elevation.csv")
    np.testing.assert_allclose(result.time, source[:, 0], rtol=0, atol=1e-12)
    np.testing.assert_allclose(result.wave_elevation, source[:, 1],
                               rtol=0, atol=1e-12)
    assert np.isfinite(result.bodies["flap"].position).all()

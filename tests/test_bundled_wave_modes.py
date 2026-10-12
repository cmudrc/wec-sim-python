"""Small wave calculations that run without a MATLAB or Applications checkout."""

from pathlib import Path

import numpy as np
import pytest
from scipy.io import savemat

from wecsim.environment.irregularWave import (
    imported_full_directional_components,
    pm_equal_energy_components,
    synthesize_full_directional_response,
    synthesize_irregular_response,
    synthesize_multiple_irregular_response,
)
from wecsim.environment.waveClass import WaveClass


ROOT = Path(__file__).resolve().parents[1]
SPHERE = ROOT / "examples/data/sphere.h5"
OSWEC = ROOT / "examples/data/oswec.h5"
WAVE_DATA = ROOT / "tests/test_objects/test_waveclass/testData"


def test_full_directional_wave_and_source_quadrature_use_distinct_force_weights(tmp_path):
    heading = np.array([0., 30., 60.])
    frequency = np.array([.1, .15, .2])
    width = np.deg2rad(30.)
    spectrum = np.array([.2, .4, .2])
    path = tmp_path / "full-directional.mat"
    savemat(path, {
        "frequencies": frequency, "spectrum": spectrum,
        "directions": heading,
        "spread": np.full((3, 3), 1 / (3 * width)),
    })
    phase = np.zeros((3, 3))
    components = imported_full_directional_components(
        OSWEC, path, phase=phase,
    )
    np.testing.assert_allclose(components.phase, phase)
    np.testing.assert_allclose(components.d_theta, width)

    integrated = synthesize_full_directional_response(
        OSWEC, components, dt=.1, end_time=.2, ramp_time=0,
    )
    matlab_weight = synthesize_full_directional_response(
        OSWEC, components, dt=.1, end_time=.2, ramp_time=0,
        force_quadrature="matlab_omitted", excitation_interpolation="spline_frequency",
    )
    energy = (components.spectral_amplitude * components.d_omega[:, None]
              * components.d_theta[None, :])
    expected = np.array([
        np.sum(np.sqrt(energy) * np.cos(t * components.omega[:, None] + phase))
        for t in integrated.time
    ])
    np.testing.assert_allclose(integrated.elevation, expected, rtol=0, atol=1e-12)
    np.testing.assert_allclose(matlab_weight.elevation, expected, rtol=0, atol=1e-12)
    assert np.isfinite(integrated.excitation_force).all()
    assert not np.allclose(integrated.excitation_force,
                           matlab_weight.excitation_force)

    savemat(path, {
        "frequencies": frequency, "spectrum": spectrum,
        "directions": heading, "spread": np.ones((3, 3)),
    })
    with pytest.raises(ValueError, match="integrate to one"):
        imported_full_directional_components(OSWEC, path, phase=phase)


def test_multiple_seas_sum_wave_and_excitation_before_dynamics():
    seas = tuple(
        pm_equal_energy_components(
            SPHERE, significant_height=height, peak_period=8,
            directions=[0], spreading=[1], count=8, seed=seed,
        )
        for height, seed in ((1.0, 1), (1.5, 2))
    )
    individual = [
        synthesize_irregular_response(
            SPHERE, sea, dt=.1, end_time=.2, ramp_time=.1,
        ) for sea in seas
    ]
    combined = synthesize_multiple_irregular_response(
        SPHERE, seas, dt=.1, end_time=.2, ramp_time=.1,
    )
    np.testing.assert_allclose(combined.elevation,
                               individual[0].elevation + individual[1].elevation)
    np.testing.assert_allclose(combined.excitation_force,
                               individual[0].excitation_force
                               + individual[1].excitation_force)


def test_legacy_imported_elevation_uses_file_values_and_ramp():
    path = WAVE_DATA / "etaImport_1_test/etaData.txt"
    source = np.loadtxt(path)
    wave = WaveClass("etaImport")
    wave.etaDataFile = str(path)
    wave.waveSetup([.05, 2], "infinite", .1, .1, 2, 9.81, 1000, .2)
    time, elevation = wave.waveAmpTime
    expected = np.interp(time, source[:, 0], source[:, 1])
    expected[0] = 0
    np.testing.assert_allclose(elevation, expected, rtol=0, atol=1e-12)
    assert np.asarray(wave.waveAmpTime1[1]).sum() == 0


def test_legacy_imported_spectrum_replays_supplied_phases(tmp_path):
    file = tmp_path / "spectrum.txt"
    frequency = np.array([.1, .15, .2])
    density = np.array([.2, .4, .2])
    phase = np.array([0., np.pi / 2, np.pi])
    np.savetxt(file, np.column_stack((frequency, density, phase)))
    wave = WaveClass("spectrumImport")
    wave.spectrumDataFile = str(file)
    wave.waveSetup([.05, 2], "infinite", 0, .1, 2, 9.81, 1000, .2)
    np.testing.assert_allclose(wave.w, 2 * np.pi * frequency)
    np.testing.assert_allclose(wave.phase[0], phase)
    expected_first = np.sum(np.sqrt(wave.A * wave.dw) * np.cos(phase))
    np.testing.assert_allclose(wave.waveAmpTime[1][0], expected_first,
                               rtol=0, atol=1e-12)

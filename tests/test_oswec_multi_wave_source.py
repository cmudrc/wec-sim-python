"""Pair each sea in the published OSWEC multiple-wave input with MATLAB."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim.environment.irregularWave import (
    pm_equal_energy_components, synthesize_irregular_response,
)


APPLICATIONS = os.environ.get("WEC_SIM_APPLICATIONS_DIR")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
pytestmark = pytest.mark.skipif(
    not (APPLICATIONS and REFERENCE
         and os.environ.get("WEC_SIM_REFERENCE_MODEL") == "OSWEC_MULTI_WAVE"),
    reason="paired MATLAB multiple-wave output not provided",
)


def test_published_settings_and_both_wave_realizations():
    reference = Path(REFERENCE)
    hydro = (Path(APPLICATIONS) / "_Common_Input_Files/OSWEC/hydroData/oswec.h5")
    settings = np.loadtxt(reference / "OSWEC_MULTI_WAVE_settings.csv", delimiter=",")
    np.testing.assert_allclose(
        settings, [12700, 1.85e6, 1.85e6, 1.85e6, 1, 0.01, 120000, 0],
        rtol=0, atol=1e-12,
    )
    components = []
    for number, height, direction in ((1, 2, 0), (2, 1, 90)):
        values = np.loadtxt(
            reference / f"OSWEC_MULTI_WAVE_wave{number}_components.csv",
            delimiter=",",
        )
        elevation = np.loadtxt(
            reference / f"OSWEC_MULTI_WAVE_wave{number}_elevation.csv",
            delimiter=",",
        )
        assert values.shape == (500, 4)
        assert elevation.shape == (1001, 2)
        assert np.isfinite(values).all() and np.isfinite(elevation).all()
        generated = pm_equal_energy_components(
            hydro, significant_height=height, peak_period=3,
            directions=np.array([direction]), spreading=np.array([1.0]),
            phase=values[:, 3, None],
        )
        np.testing.assert_allclose(generated.omega, values[:, 0], rtol=0, atol=1e-12)
        np.testing.assert_allclose(generated.spectral_amplitude, values[:, 1],
                                   rtol=0, atol=1e-12)
        np.testing.assert_allclose(generated.d_omega, values[:, 2], rtol=0, atol=1e-12)
        response = synthesize_irregular_response(
            hydro, generated, dt=0.1, end_time=100, ramp_time=0.1,
        )
        np.testing.assert_allclose(response.time, elevation[:, 0], rtol=0, atol=1e-10)
        assert np.max(np.abs(response.elevation - elevation[:, 1])) < 1e-11
        components.append(values)

    np.testing.assert_allclose(components[0][:, 0], components[1][:, 0],
                               rtol=0, atol=1e-12)
    np.testing.assert_allclose(components[0][:, 1], 4 * components[1][:, 1],
                               rtol=0, atol=1e-12)
    for number in (1, 2):
        body = np.loadtxt(
            reference / f"OSWEC_MULTI_WAVE_Multiple_Wave_Spectra_body{number}.csv",
            delimiter=",",
        )
        assert body.shape == (1001, 25)
        assert np.isfinite(body).all()
        if number == 1:
            # The published Simscape hinge pitches about y. Passive-yaw
            # interpolation is enabled, but the body has no yaw motion.
            assert np.max(np.abs(body[:, 6])) < 1e-10
        else:
            assert np.max(np.abs(body[:, 1:7] - body[0, 1:7])) < 1e-10
    pto = np.loadtxt(
        reference / "OSWEC_MULTI_WAVE_Multiple_Wave_Spectra_pto1.csv",
        delimiter=",",
    )
    assert pto.shape == (1001, 25)
    assert np.isfinite(pto).all()

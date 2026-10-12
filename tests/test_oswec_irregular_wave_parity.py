"""Rebuild the current OSWEC wave and excitation from paired MATLAB phases."""

import os
from pathlib import Path
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from wecsim.dynamics.hingePitch import solve_hinged_pitch_from_excitation  # noqa: E402
from wecsim.environment.irregularWave import (  # noqa: E402
    IrregularComponents, pm_equal_energy_components,
    synthesize_irregular_response,
)

CORE = os.environ.get("WEC_SIM_MATLAB_CORE_DIR")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
pytestmark = pytest.mark.skipif(
    not (CORE and REFERENCE and os.environ.get("WEC_SIM_REFERENCE_MODEL") == "OSWEC"),
    reason="paired MATLAB OSWEC wave components not provided",
)


@pytest.fixture(scope="module")
def comparison():
    reference = Path(REFERENCE)
    hydro_file = Path(CORE) / "examples" / "OSWEC" / "hydroData" / "oswec.h5"
    values = np.loadtxt(reference / "OSWEC_wave_components.csv", delimiter=",")
    directions = np.loadtxt(reference / "OSWEC_wave_directions.csv", delimiter=",")
    matlab_wave = np.loadtxt(reference / "OSWEC_wave_elevation.csv", delimiter=",")
    matlab_body = np.loadtxt(reference / "OSWEC_OSWEC_body1.csv", delimiter=",")
    matlab_components = IrregularComponents(
        omega=values[:, 0], spectral_amplitude=values[:, 1],
        d_omega=values[:, 2], directions=directions[:, 0],
        spreading=directions[:, 1], phase=values[:, 3:],
    )
    dt = matlab_wave[1, 0] - matlab_wave[0, 0]
    generated = pm_equal_energy_components(
        hydro_file, significant_height=2.5, peak_period=8,
        directions=matlab_components.directions,
        spreading=matlab_components.spreading, phase=matlab_components.phase,
    )
    response = synthesize_irregular_response(
        hydro_file, generated, dt=dt, end_time=matlab_wave[-1, 0],
        ramp_time=100,
    )
    return hydro_file, matlab_components, generated, response, matlab_wave, matlab_body


def test_oswec_pm_equal_energy_bins(comparison):
    _, components, generated, _, _, _ = comparison
    np.testing.assert_allclose(generated.omega, components.omega, rtol=0, atol=1e-12)
    np.testing.assert_allclose(generated.d_omega, components.d_omega, rtol=0, atol=1e-12)
    np.testing.assert_allclose(
        generated.spectral_amplitude, components.spectral_amplitude, rtol=0, atol=1e-12,
    )


def test_oswec_wave_elevation_and_excitation(comparison):
    _, _, _, response, matlab_wave, matlab_body = comparison
    np.testing.assert_allclose(response.time, matlab_wave[:, 0], rtol=0, atol=1e-10)
    assert np.max(np.abs(response.elevation - matlab_wave[:, 1])) < 1e-11
    assert np.max(np.abs(response.excitation_force - matlab_body[:, 19:25])) < 1e-6


def test_oswec_python_excitation_drives_matlab_matched_pitch(comparison):
    hydro_file, _, _, wave, _, matlab_body = comparison
    motion = solve_hinged_pitch_from_excitation(
        hydro_file, wave.excitation_force, hinge_z=-8.9,
        body_mass=127_000, pitch_inertia=1.85e6, pto_damping=12_000,
    )
    assert np.max(np.abs(motion.angle - matlab_body[:, 5])) < 0.004

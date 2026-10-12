"""Check the OSWEC hinged-pitch subsystem against paired MATLAB output."""

import os
from pathlib import Path
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from wecsim.dynamics.hingePitch import solve_hinged_pitch_from_excitation  # noqa: E402

CORE = os.environ.get("WEC_SIM_MATLAB_CORE_DIR")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
pytestmark = pytest.mark.skipif(
    not (CORE and REFERENCE and os.environ.get("WEC_SIM_REFERENCE_MODEL") == "OSWEC"),
    reason="paired MATLAB OSWEC excitation and response not provided",
)


@pytest.fixture(scope="module")
def comparison():
    expected = np.loadtxt(
        Path(REFERENCE) / "OSWEC_OSWEC_body1.csv", delimiter=",",
    )
    h5_file = Path(CORE) / "examples" / "OSWEC" / "hydroData" / "oswec.h5"
    response = solve_hinged_pitch_from_excitation(
        h5_file, expected[:, 19:25], hinge_z=-8.9,
        body_mass=127_000, pitch_inertia=1.85e6, pto_damping=12_000,
    )
    return response, expected


def test_oswec_hinged_pitch_against_matlab(comparison):
    response, expected = comparison
    np.testing.assert_allclose(response.time, expected[:, 0], rtol=0, atol=1e-10)
    assert np.max(np.abs(response.angle - expected[:, 5])) < 0.004
    assert np.max(np.abs(response.angular_velocity - expected[:, 11])) < 0.005
    assert np.max(np.abs(response.center_position[:, 0] - expected[:, 1])) < 0.02
    assert np.max(np.abs(response.center_position[:, 2] - expected[:, 3])) < 0.005


def test_oswec_pto_against_matlab(comparison):
    response, _ = comparison
    expected = np.loadtxt(
        Path(REFERENCE) / "OSWEC_OSWEC_pto1.csv", delimiter=",",
    )
    np.testing.assert_allclose(response.time, expected[:, 0], rtol=0, atol=1e-10)
    assert np.max(np.abs(response.angle - expected[:, 5])) < 0.004
    assert np.max(np.abs(response.angular_velocity - expected[:, 11])) < 0.005
    assert np.max(np.abs(response.pto_torque - expected[:, 17])) < 60

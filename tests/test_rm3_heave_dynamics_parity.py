"""Compare RM3's heave subsystem with the pinned MATLAB reference run."""

import os
from pathlib import Path
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from wecsim.dynamics.linearHeave import solve_two_body_regular_heave  # noqa: E402

CORE = os.environ.get("WEC_SIM_MATLAB_CORE_DIR")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
pytestmark = pytest.mark.skipif(
    not (CORE and REFERENCE and os.environ.get("WEC_SIM_REFERENCE_MODEL") == "RM3"),
    reason="live MATLAB RM3 inputs and body trajectories not provided",
)


@pytest.fixture(scope="module")
def response():
    h5_file = Path(CORE) / "examples" / "RM3" / "hydroData" / "rm3.h5"
    return solve_two_body_regular_heave(
        h5_file, wave_height=2.5, wave_period=8.0, pto_damping=1_200_000.0,
    )


@pytest.mark.parametrize("body_number", [1, 2])
def test_rm3_regular_wave_heave_against_matlab(response, body_number):
    expected = np.loadtxt(
        Path(REFERENCE) / f"RM3_RM3_body{body_number}.csv", delimiter=",",
    )
    index = body_number - 1
    np.testing.assert_allclose(response.time, expected[:, 0], rtol=0, atol=1e-10)
    # Body-position columns 1:7; body-velocity columns 7:13. Heave is index 2.
    assert np.max(np.abs(response.position[:, index] - expected[:, 3])) < 0.0065
    assert np.max(np.abs(response.velocity[:, index] - expected[:, 9])) < 0.0065
    # Excitation is an exact use of the same hydrodynamic coefficients and
    # half-cosine ramp. This checks its phase and units separately from motion.
    np.testing.assert_allclose(
        response.excitation_force[:, index], expected[:, 21], rtol=0, atol=1e-5,
    )


def test_rm3_heave_pto_against_matlab(response):
    expected = np.loadtxt(
        Path(REFERENCE) / "RM3_RM3_pto1.csv", delimiter=",",
    )
    relative_displacement = response.displacement[:, 0] - response.displacement[:, 1]
    relative_velocity = response.velocity[:, 0] - response.velocity[:, 1]
    np.testing.assert_allclose(response.time, expected[:, 0], rtol=0, atol=1e-10)
    assert np.max(np.abs(relative_displacement - expected[:, 3])) < 0.009
    assert np.max(np.abs(relative_velocity - expected[:, 9])) < 0.0075
    assert np.max(np.abs(response.pto_force - expected[:, 15])) < 9_000
    assert np.max(np.abs(response.pto_force * relative_velocity - expected[:, 21])) < 12_000

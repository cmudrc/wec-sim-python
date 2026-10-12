"""Compare the coupled RM3 reference model's active DOFs and PTO with MATLAB."""

import os
from pathlib import Path
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from wecsim.reference_models.rm3Regular import solve_rm3_regular  # noqa: E402

CORE = os.environ.get("WEC_SIM_MATLAB_CORE_DIR")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
pytestmark = pytest.mark.skipif(
    not (CORE and REFERENCE and os.environ.get("WEC_SIM_REFERENCE_MODEL") == "RM3"),
    reason="current MATLAB RM3 input and response not provided",
)


@pytest.fixture(scope="module")
def response():
    hydro = Path(CORE) / "examples" / "RM3" / "hydroData" / "rm3.h5"
    return solve_rm3_regular(hydro)


@pytest.mark.parametrize("body_number", [1, 2])
def test_rm3_surge_heave_pitch_against_matlab(response, body_number):
    expected = np.loadtxt(
        Path(REFERENCE) / f"RM3_RM3_body{body_number}.csv", delimiter=",",
    )
    position = response.body_position[:, body_number - 1, :]
    velocity = response.body_velocity[:, body_number - 1, :]
    np.testing.assert_allclose(response.time, expected[:, 0], rtol=0, atol=1e-10)
    for dof, tolerance in ((0, 0.04), (2, 0.005), (4, 0.001)):
        assert np.max(np.abs(position[:, dof] - expected[:, 1 + dof])) < tolerance
    for dof, tolerance in ((0, 0.035), (2, 0.004), (4, 0.001)):
        assert np.max(np.abs(velocity[:, dof] - expected[:, 7 + dof])) < tolerance


def test_rm3_coupled_pto_force_against_matlab(response):
    expected = np.loadtxt(
        Path(REFERENCE) / "RM3_RM3_pto1.csv", delimiter=",",
    )
    np.testing.assert_allclose(response.time, expected[:, 0], rtol=0, atol=1e-10)
    assert np.max(np.abs(response.pto_force - expected[:, 15])) < 5_000

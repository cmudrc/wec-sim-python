"""Compare both published RM3 body-to-body application cases with MATLAB."""

import os
from pathlib import Path
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from wecsim.dynamics.rm3Regular import solve_rm3_regular  # noqa: E402

APPLICATIONS = os.environ.get("WEC_SIM_APPLICATIONS_DIR")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
pytestmark = pytest.mark.skipif(
    not (APPLICATIONS and REFERENCE),
    reason="paired MATLAB RM3 body-to-body application output not provided",
)


@pytest.mark.parametrize("case,b2b", [("B2B_Case1", False), ("B2B_Case2", True)])
def test_rm3_body_to_body_cases(case, b2b):
    hydro = Path(APPLICATIONS) / "_Common_Input_Files" / "RM3" / "hydroData" / "rm3.h5"
    response = solve_rm3_regular(hydro, b2b=b2b)
    for body_number in (1, 2):
        expected = np.loadtxt(
            Path(REFERENCE) / f"RM3_B2B_{case}_body{body_number}.csv", delimiter=",",
        )
        np.testing.assert_allclose(response.time, expected[:, 0], rtol=0, atol=1e-10)
        position = response.body_position[:, body_number - 1, :]
        velocity = response.body_velocity[:, body_number - 1, :]
        for dof, tolerance in ((0, 0.04), (2, 0.005), (4, 0.0012)):
            assert np.max(np.abs(position[:, dof] - expected[:, 1 + dof])) < tolerance
        for dof, tolerance in ((0, 0.035), (2, 0.004), (4, 0.001)):
            assert np.max(np.abs(velocity[:, dof] - expected[:, 7 + dof])) < tolerance
    expected_pto = np.loadtxt(
        Path(REFERENCE) / f"RM3_B2B_{case}_pto1.csv", delimiter=",",
    )
    np.testing.assert_allclose(response.time, expected_pto[:, 0], rtol=0, atol=1e-10)
    assert np.max(np.abs(response.pto_force - expected_pto[:, 15])) < 6_000
    if b2b:
        uncoupled = solve_rm3_regular(hydro, b2b=False)
        expected_heave = np.loadtxt(
            Path(REFERENCE) / f"RM3_B2B_{case}_body1.csv", delimiter=",",
        )[:, 3]
        coupled_error = np.sqrt(np.mean(
            (response.body_position[:, 0, 2] - expected_heave)**2,
        ))
        uncoupled_error = np.sqrt(np.mean(
            (uncoupled.body_position[:, 0, 2] - expected_heave)**2,
        ))
        assert coupled_error < uncoupled_error / 4

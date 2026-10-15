"""Compare both published RM3 body-to-body application cases with MATLAB."""

import os
from pathlib import Path
import sys

import numpy as np

import pytest
from scipy.integrate import trapezoid

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from wecsim import RegularWave, WEC  # noqa: E402
from wecsim.dynamics.rm3Regular import solve_rm3_regular  # noqa: E402

APPLICATIONS = os.environ.get("WEC_SIM_APPLICATIONS_DIR")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
PORTABLE = ROOT / "tests/fixtures/matlab/rm3_b2b_cic/source.npz"
HYDRO = ROOT / "examples/data/rm3.h5"


def _source_array(case, signal):
    if APPLICATIONS and REFERENCE:
        return np.loadtxt(
            Path(REFERENCE) / f"RM3_B2B_{case}_{signal}.csv", delimiter=",",
        )
    with np.load(PORTABLE, allow_pickle=False) as saved:
        return saved[f"{case}_{signal}"].copy()


@pytest.mark.parametrize("case,b2b", [("B2B_Case1", False), ("B2B_Case2", True)])
def test_rm3_body_to_body_cases(case, b2b):
    hydro = (Path(APPLICATIONS) / "_Common_Input_Files/RM3/hydroData/rm3.h5"
             if APPLICATIONS and REFERENCE else HYDRO)
    response = solve_rm3_regular(hydro, b2b=b2b)
    for body_number in (1, 2):
        expected = _source_array(case, f"body{body_number}")
        np.testing.assert_allclose(response.time, expected[:, 0], rtol=0, atol=1e-10)
        position = response.body_position[:, body_number - 1, :]
        velocity = response.body_velocity[:, body_number - 1, :]
        for dof, tolerance in ((0, 0.04), (2, 0.005), (4, 0.0012)):
            assert np.max(np.abs(position[:, dof] - expected[:, 1 + dof])) < tolerance
        for dof, tolerance in ((0, 0.035), (2, 0.004), (4, 0.001)):
            assert np.max(np.abs(velocity[:, dof] - expected[:, 7 + dof])) < tolerance
    expected_pto = _source_array(case, "pto1")
    np.testing.assert_allclose(response.time, expected_pto[:, 0], rtol=0, atol=1e-10)
    assert np.max(np.abs(response.pto_force - expected_pto[:, 15])) < 6_000
    if b2b:
        uncoupled = solve_rm3_regular(hydro, b2b=False)
        expected_heave = _source_array(case, "body1")[:, 3]
        coupled_error = np.sqrt(np.mean(
            (response.body_position[:, 0, 2] - expected_heave)**2,
        ))
        uncoupled_error = np.sqrt(np.mean(
            (uncoupled.body_position[:, 0, 2] - expected_heave)**2,
        ))
        assert coupled_error < uncoupled_error / 4


@pytest.mark.parametrize("case,b2b", [("B2B_Case1", False),
                                      ("B2B_Case2", True)])
def test_public_rm3_regular_body_to_body_against_matlab(case, b2b):
    """Gate the installable WEC builder against both published regular cases."""
    hydro = (Path(APPLICATIONS) / "_Common_Input_Files/RM3/hydroData/rm3.h5"
             if APPLICATIONS and REFERENCE else HYDRO)
    wec = WEC(case, body_to_body=b2b)
    float_body = wec.body("float", hydro, inertia=(0, 21_306_090.66, 0))
    spar = wec.body("spar", hydro, inertia=(0, 94_407_091.24, 0))
    wec.floating_joint(float_body, spar, damping=1_200_000)
    response = wec.run(RegularWave(2.5, 8), dt=.1, end_time=400,
                       ramp_time=100)

    for body_number, name in ((1, "float"), (2, "spar")):
        expected = _source_array(case, f"body{body_number}")
        np.testing.assert_allclose(response.time, expected[:, 0], rtol=0,
                                   atol=1e-10)
        body = response.bodies[name]
        for dof, position_limit, speed_limit in (
            (0, .04, .035), (2, .005, .005), (4, .0012, .001),
        ):
            assert np.max(np.abs(body.position[:, dof]
                                 - expected[:, 1 + dof])) < position_limit
            assert np.max(np.abs(body.velocity[:, dof]
                                 - expected[:, 7 + dof])) < speed_limit

    expected_pto = _source_array(case, "pto1")
    pto = response.ptos["relative_heave"]
    np.testing.assert_allclose(response.time, expected_pto[:, 0], rtol=0,
                               atol=1e-10)
    assert np.max(np.abs(pto.stroke - expected_pto[:, 3])) < .008
    assert np.max(np.abs(pto.velocity - expected_pto[:, 9])) < .006
    assert np.max(np.abs(pto.force - expected_pto[:, 15])) < 6_000
    # MATLAB reports mechanical PTO output with the opposite sign to Python's
    # absorbed-power channel. Integrate both to catch a persistent bias.
    assert np.max(np.abs(pto.absorbed_power + expected_pto[:, 21])) < 8_000
    source_work = -trapezoid(expected_pto[:, 21], x=response.time)
    python_work = trapezoid(pto.absorbed_power, x=response.time)
    assert abs(python_work - source_work) / source_work < .01

"""Compare published RM3 regularCIC body-to-body cases with MATLAB."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim import RegularCICWave, WEC
from wecsim.dynamics.caseDynamics import run_case
from wecsim.dynamics.rm3Regular import solve_rm3_regular


APPLICATIONS = os.environ.get("WEC_SIM_APPLICATIONS_DIR")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
PORTABLE = Path(__file__).parent / "fixtures/matlab/rm3_b2b_cic/source.npz"
HYDRO = Path(__file__).resolve().parents[1] / "examples/data/rm3.h5"
external_only = pytest.mark.skipif(
    not (APPLICATIONS and REFERENCE),
    reason="paired MATLAB RM3 regularCIC output not provided",
)


def _source_array(case, signal):
    if APPLICATIONS and REFERENCE:
        return np.loadtxt(
            Path(REFERENCE) / f"RM3_B2B_{case}_{signal}.csv", delimiter=",",
        )
    with np.load(PORTABLE, allow_pickle=False) as saved:
        return saved[f"{case}_{signal}"].copy()


def _max_error(actual, expected, limit, label):
    actual = np.asarray(actual)
    expected = np.asarray(expected)
    assert actual.shape == expected.shape, label
    assert np.isfinite(expected).all(), label
    error = np.max(np.abs(actual - expected))
    assert error < limit, f"{label}: max error {error:.6g} exceeds {limit}"


@pytest.mark.parametrize("case,b2b", [
    ("B2B_Case3", False), ("B2B_Case4", True),
])
@pytest.mark.parametrize("mass_scheme", [None, "simulink_delay"],
                         ids=["default-implicit", "source-delay"])
def test_rm3_regular_cic_against_matlab(case, b2b, mass_scheme):
    """Gate the ordinary solver as well as the opt-in Simulink numerics."""
    hydro = ((Path(APPLICATIONS) / "_Common_Input_Files/RM3/hydroData/rm3.h5").resolve()
             if APPLICATIONS and REFERENCE else HYDRO)
    simulation = {"dt": 0.1, "end_time": 400, "ramp_time": 100,
                  "radiation_memory": 60}
    if mass_scheme is not None:
        simulation["added_mass_scheme"] = mass_scheme
    config = {
        "simulation": simulation,
        "wave": {"type": "regularCIC", "height": 2.5, "period": 8},
        "bodies": [
            {"hydro_file": str(hydro), "hydro_body": 1,
             "mass": "equilibrium", "pitch_inertia": 21_306_090.66},
            {"hydro_file": str(hydro), "hydro_body": 2,
             "mass": "equilibrium", "pitch_inertia": 94_407_091.24},
        ],
        "constraint": {"kind": "floating_joint", "location": [0, 0, 0]},
        "pto": {"kind": "relative_heave", "damping": 1_200_000},
        "body_to_body": b2b,
    }
    response = run_case(config)
    for body in (1, 2):
        expected = _source_array(case, f"body{body}")
        np.testing.assert_allclose(response.time, expected[:, 0], rtol=0, atol=1e-10)
        for dof, label, position_limit, velocity_limit in (
            (0, "surge", 0.075, 0.04),
            (2, "heave", 0.005, 0.0035),
            (4, "pitch", 0.002 if mass_scheme is None else 0.0015,
             0.0011),
        ):
            _max_error(response.body_position[:, body - 1, dof],
                       expected[:, 1 + dof], position_limit,
                       f"{case} body{body} {label} position")
            _max_error(response.body_velocity[:, body - 1, dof],
                       expected[:, 7 + dof], velocity_limit,
                       f"{case} body{body} {label} velocity")
    pto = _source_array(case, "pto1")
    np.testing.assert_allclose(response.time, pto[:, 0], rtol=0, atol=1e-10)
    # Recover the PTO translation at the floating joint, accounting for
    # the two body-center offsets when the assembly pitches.
    center_gap = (response.body_position[0, 0, 2]
                  - response.body_position[0, 1, 2])
    pitch = response.body_position[:, 0, 4]
    pitch_speed = response.body_velocity[:, 0, 4]
    cosine = np.cos(pitch)
    stroke = ((response.body_position[:, 0, 2]
               - response.body_position[:, 1, 2]) / cosine - center_gap)
    speed = ((response.body_velocity[:, 0, 2]
              - response.body_velocity[:, 1, 2]
              + (center_gap + stroke) * np.sin(pitch) * pitch_speed) / cosine)
    _max_error(stroke, pto[:, 3], 0.008, f"{case} PTO stroke")
    _max_error(speed, pto[:, 9], 0.005, f"{case} PTO speed")
    _max_error(response.pto_force, pto[:, 15], 6_000, f"{case} PTO force")
    _max_error(response.pto_force * speed, pto[:, 21], 7_000,
               f"{case} PTO power")


@external_only
@pytest.mark.parametrize("convolution_case,fitted_case,b2b,source_gap", [
    ("B2B_Case3", "B2B_Case5", False, (0.39, 0.42)),
    ("B2B_Case4", "B2B_Case6", True, (0.59, 0.63)),
])
def test_state_space_source_motion_differs_from_convolution(
    convolution_case, fitted_case, b2b, source_gap,
):
    """Keep the physical convolution comparison separate from a nonpassive fit.

    The published case pair changes only its radiation option. This gate
    compares complete source trajectories as well as Python's unchanged
    physical convolution path; it does not claim fitted-case parity.
    """
    hydro = (Path(APPLICATIONS) / "_Common_Input_Files/RM3/hydroData/rm3.h5")
    reference = Path(REFERENCE)
    convolution = np.loadtxt(
        reference / f"RM3_B2B_{convolution_case}_body1.csv", delimiter=",",
    )
    fitted = np.loadtxt(
        reference / f"RM3_B2B_{fitted_case}_body1.csv", delimiter=",",
    )
    assert convolution.shape == fitted.shape == (4001, 25)
    np.testing.assert_allclose(convolution[:, 0], fitted[:, 0], rtol=0,
                               atol=1e-10)
    physical = solve_rm3_regular(hydro, b2b=b2b, radiation_memory=60)
    np.testing.assert_allclose(physical.time, convolution[:, 0], rtol=0,
                               atol=1e-10)
    source_difference = np.max(np.abs(fitted[:, 1] - convolution[:, 1]))
    assert source_gap[0] < source_difference < source_gap[1]
    physical_difference = np.max(np.abs(
        physical.body_position[:, 0, 0] - convolution[:, 1],
    ))
    assert physical_difference < 0.075
    fitted_difference = np.max(np.abs(
        physical.body_position[:, 0, 0] - fitted[:, 1],
    ))
    assert fitted_difference > 0.35


@external_only
@pytest.mark.parametrize("case,b2b", [
    ("B2B_Case5", False), ("B2B_Case6", True),
])
def test_opt_in_fitted_radiation_advances_published_motion(case, b2b):
    """Pair the nonpassive source fit without changing convolution dynamics."""
    hydro = (Path(APPLICATIONS) / "_Common_Input_Files/RM3/hydroData/rm3.h5")
    wec = WEC(f"RM3 {case} source fit", body_to_body=b2b)
    float_body = wec.body("float", hydro, inertia=(0, 21_306_090.66, 0))
    spar = wec.body("spar", hydro, inertia=(0, 94_407_091.24, 0))
    wec.floating_joint(float_body, spar, damping=1_200_000,
                       radiation_method="state_space")
    response = wec.run(
        RegularCICWave(2.5, 8), dt=.1, end_time=400,
        ramp_time=100, radiation_memory=60,
    )
    assert response.case["simulation"]["radiation_method"] == "state_space"
    reference = Path(REFERENCE)
    for body_index, name in enumerate(("float", "spar"), start=1):
        expected = np.loadtxt(
            reference / f"RM3_B2B_{case}_body{body_index}.csv",
            delimiter=",",
        )
        _max_error(response.time, expected[:, 0], 1e-10, "time")
        for dof, label, position_limit, velocity_limit in (
            (0, "surge", .08, .04),
            (2, "heave", .004, .0035),
            (4, "pitch", .0025, .0015),
        ):
            _max_error(response.bodies[name].position[:, dof],
                       expected[:, 1 + dof], position_limit,
                       f"{case} {name} {label} position")
            _max_error(response.bodies[name].velocity[:, dof],
                       expected[:, 7 + dof], velocity_limit,
                       f"{case} {name} {label} velocity")
    expected_pto = np.loadtxt(reference / f"RM3_B2B_{case}_pto1.csv",
                              delimiter=",")
    pto = response.ptos["relative_heave"]
    _max_error(pto.stroke, expected_pto[:, 3], .006, f"{case} PTO stroke")
    _max_error(pto.velocity, expected_pto[:, 9], .004, f"{case} PTO speed")
    _max_error(pto.force, expected_pto[:, 15], 5_000, f"{case} PTO force")
    _max_error(pto.absorbed_power, -expected_pto[:, 21], 6_000,
               f"{case} absorbed power")

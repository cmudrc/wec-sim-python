"""Pair the published OSWEC nonhydrodynamic-base case with MATLAB output."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim import RegularWave, WEC, WorldPoint
from wecsim.dynamics.caseDynamics import run_case
from wecsim.dynamics.hingePitch import solve_hinged_pitch_regular


APPLICATIONS = os.environ.get("WEC_SIM_APPLICATIONS_DIR")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
pytestmark = pytest.mark.skipif(
    not (APPLICATIONS and REFERENCE),
    reason="paired MATLAB OSWEC nonhydrodynamic output not provided",
)


def _max_error(actual, expected, limit, label):
    assert actual.shape == expected.shape, label
    assert np.isfinite(expected).all(), label
    error = np.max(np.abs(actual - expected))
    assert error < limit, f"{label}: max error {error:.6g} exceeds {limit}"


def test_published_nonhydrodynamic_base_against_matlab():
    hydro = (Path(APPLICATIONS)
             / "_Common_Input_Files/OSWEC/hydroData/oswec.h5").resolve()
    case = {
        "simulation": {"dt": 0.1, "end_time": 400, "ramp_time": 100},
        "wave": {"type": "regular", "height": 2.5, "period": 8},
        "bodies": [
            {"name": "flap", "hydro_file": str(hydro), "hydro_body": 1,
             "mass": 127_000, "pitch_inertia": 1.85e6},
            {"name": "base", "nonhydro": True, "fixed": True,
             "center_gravity": [0, 0, -10.9], "mass": 999,
             "inertia": [1, 1, 1]},
        ],
        "constraint": {"kind": "fixed_hinge", "location": [0, 0, -10]},
        "pto": {"kind": "pitch", "damping": 0,
                "location": [0, 0, -8.9]},
    }
    response = run_case(case)
    reference = Path(REFERENCE)
    flap = np.loadtxt(
        reference / "OSWEC_Nonhydro_Nonhydro_body1.csv", delimiter=",",
    )
    base = np.loadtxt(
        reference / "OSWEC_Nonhydro_Nonhydro_body2.csv", delimiter=",",
    )
    np.testing.assert_allclose(response.time, flap[:, 0], rtol=0, atol=1e-10)
    np.testing.assert_allclose(response.time, base[:, 0], rtol=0, atol=1e-10)
    for dof, label, position_limit, velocity_limit in (
        (0, "surge", 0.03, 0.025),
        (2, "heave", 0.012, 0.01),
        (4, "pitch", 0.007, 0.006),
    ):
        _max_error(response.body_position[:, 0, dof], flap[:, 1 + dof],
                   position_limit, f"flap {label} position")
        _max_error(response.body_velocity[:, 0, dof], flap[:, 7 + dof],
                   velocity_limit, f"flap {label} velocity")
    for dof in (1, 3, 5):
        _max_error(response.body_position[:, 0, dof], flap[:, 1 + dof],
                   1e-10, f"flap stationary DOF {dof} position")
        _max_error(response.body_velocity[:, 0, dof], flap[:, 7 + dof],
                   1e-10, f"flap stationary DOF {dof} velocity")
    np.testing.assert_allclose(
        response.body_position[:, 1, :], base[:, 1:7], rtol=0, atol=1e-10,
    )
    np.testing.assert_allclose(
        response.body_velocity[:, 1, :], base[:, 7:13], rtol=0, atol=1e-10,
    )
    pto = np.loadtxt(
        reference / "OSWEC_Nonhydro_Nonhydro_pto1.csv", delimiter=",",
    )
    _max_error(response.body_position[:, 0, 4], pto[:, 5], 0.007, "PTO angle")
    _max_error(response.body_velocity[:, 0, 4], pto[:, 11], 0.006,
               "PTO angular speed")
    _max_error(response.pto_force, pto[:, 17], 1e-6, "PTO torque")
    direct = solve_hinged_pitch_regular(
        hydro, wave_height=2.5, wave_period=8, hinge_z=-8.9,
        body_mass=127_000, pitch_inertia=1.85e6, pto_damping=0,
    )
    _max_error(direct.excitation_force, flap[:, 19:25], 1.0,
               "six-component flap excitation force")

    configured = WEC("OSWEC nonhydrodynamic base")
    moving = configured.body("flap", hydro, mass=127_000,
                             inertia=(0, 1.85e6, 0))
    fixed = configured.fixed_body(
        "base", center_gravity=(0, 0, -10.9), mass=999,
        inertia=(1, 1, 1),
    )
    configured.fixed_hinge(
        moving, fixed, location=WorldPoint(0, 0, -10),
        pto_location=WorldPoint(0, 0, -8.9), damping=0,
    )
    public = configured.run(
        RegularWave(2.5, 8), dt=0.1, end_time=400, ramp_time=100,
    )
    np.testing.assert_allclose(public.raw.body_position,
                               response.body_position, rtol=0, atol=1e-12)
    np.testing.assert_allclose(public.raw.body_velocity,
                               response.body_velocity, rtol=0, atol=1e-12)
    _max_error(public.bodies["flap"].position[:, 4], flap[:, 5],
               0.007, "configured flap pitch")
    np.testing.assert_allclose(public.bodies["base"].position,
                               base[:, 1:7], rtol=0, atol=1e-10)
    np.testing.assert_allclose(public.ptos["hinge"].force,
                               pto[:, 17], rtol=0, atol=1e-6)
    np.testing.assert_allclose(public.coordinates["pitch"].position,
                               public.ptos["hinge"].stroke, rtol=0, atol=0)

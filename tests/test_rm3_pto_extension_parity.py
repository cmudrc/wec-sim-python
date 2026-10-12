"""Compare both published RM3 PTO-extension free decays with MATLAB."""

import os
from pathlib import Path

import h5py
import numpy as np
import pytest

from wecsim.dynamics.caseDynamics import run_case


APPLICATIONS = os.environ.get("WEC_SIM_APPLICATIONS_DIR")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
pytestmark = pytest.mark.skipif(
    not (APPLICATIONS and REFERENCE),
    reason="paired MATLAB RM3 PTO-extension output not provided",
)


def _max_error(actual, expected, limit, label):
    assert actual.shape == expected.shape, label
    assert np.isfinite(expected).all(), label
    error = np.max(np.abs(actual - expected))
    assert error < limit, f"{label}: max error {error:.6g} exceeds {limit}"


@pytest.mark.parametrize("name,coordinate,active_body,heave_limit,speed_limit", [
    ("float", {"float_heave": 5}, 0, 0.090, 0.210),
    ("spar", {"spar_heave": -5}, 1, 0.018, 0.008),
])
def test_published_pto_extension_free_decay(
    name, coordinate, active_body, heave_limit, speed_limit,
):
    hydro = (Path(APPLICATIONS)
             / "_Common_Input_Files/RM3/hydroData/rm3.h5").resolve()
    case = {
        "simulation": {"dt": 0.1, "end_time": 30,
                       "radiation_memory": 60,
                       "added_mass_scheme": "simulink_delay"},
        "wave": {"type": "none"},
        "bodies": [
            {"hydro_file": str(hydro), "hydro_body": 1,
             "mass": "equilibrium", "pitch_inertia": 21_306_090.66},
            {"hydro_file": str(hydro), "hydro_body": 2,
             "mass": "equilibrium", "pitch_inertia": 94_407_091.24},
        ],
        "constraint": {"kind": "floating_joint", "location": [0, 0, 0],
                       "initial_coordinate": coordinate},
        "pto": {"kind": "relative_heave", "damping": 0},
    }
    response = run_case(case)
    reference = Path(REFERENCE)
    assert response.wave_elevation is None
    assert response.body_position.shape == (301, 2, 6)
    for index in range(2):
        expected = np.loadtxt(
            reference / f"RM3_PTO_Extension_{name}_body{index+1}.csv",
            delimiter=",",
        )
        np.testing.assert_allclose(response.time, expected[:, 0],
                                   rtol=0, atol=1e-10)
        for dof, axis, position_limit, velocity_limit in (
            (0, "surge", 3e-5, 1e-5),
            (2, "heave", heave_limit if index == active_body else 1e-8,
             speed_limit if index == active_body else 1e-8),
            (4, "pitch", 1e-6, 1e-6),
        ):
            _max_error(response.body_position[:, index, dof],
                       expected[:, 1 + dof], position_limit,
                       f"{name} body{index+1} {axis} position")
            _max_error(response.body_velocity[:, index, dof],
                       expected[:, 7 + dof], velocity_limit,
                       f"{name} body{index+1} {axis} velocity")
        for dof in (1, 3, 5):
            _max_error(response.body_position[:, index, dof],
                       expected[:, 1 + dof], 1e-9,
                       f"{name} body{index+1} stationary position {dof}")
            _max_error(response.body_velocity[:, index, dof],
                       expected[:, 7 + dof], 1e-9,
                       f"{name} body{index+1} stationary velocity {dof}")

    with h5py.File(hydro) as file:
        center_gap = (
            np.asarray(file["body1/properties/cg"]).ravel()[2]
            - np.asarray(file["body2/properties/cg"]).ravel()[2]
        )
    pitch = response.body_position[:, 0, 4]
    pitch_speed = response.body_velocity[:, 0, 4]
    cosine = np.cos(pitch)
    stroke = ((response.body_position[:, 0, 2]
               - response.body_position[:, 1, 2]) / cosine - center_gap)
    speed = ((response.body_velocity[:, 0, 2]
              - response.body_velocity[:, 1, 2]
              + (center_gap + stroke) * np.sin(pitch) * pitch_speed) / cosine)
    pto = np.loadtxt(
        reference / f"RM3_PTO_Extension_{name}_pto1.csv", delimiter=",",
    )
    np.testing.assert_allclose(response.time, pto[:, 0], rtol=0, atol=1e-10)
    _max_error(stroke, pto[:, 3], heave_limit, f"{name} PTO stroke")
    _max_error(speed, pto[:, 9], speed_limit, f"{name} PTO speed")
    _max_error(response.pto_force, pto[:, 15], 1e-7,
               f"{name} zero PTO force")
    _max_error(response.pto_force * speed, pto[:, 21], 1e-7,
               f"{name} zero PTO power")

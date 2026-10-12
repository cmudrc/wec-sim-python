"""Pair RM3 PTO-extension free decay at three pinned MATLAB solver steps."""

import os
from pathlib import Path

import h5py
import numpy as np
import pytest

from wecsim.dynamics.caseDynamics import run_case


APPLICATIONS = os.environ.get("WEC_SIM_APPLICATIONS_DIR")
SOURCE = os.environ.get("WEC_SIM_MATLAB_STEP_AUDIT_DIR")
pytestmark = pytest.mark.skipif(
    not (APPLICATIONS and SOURCE),
    reason="pinned RM3 PTO-extension step audit not provided",
)


def _python_response(hydro, coordinate, dt):
    return run_case({
        "simulation": {"dt": dt, "end_time": 30,
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
    })


@pytest.mark.parametrize("name,coordinate,active,position_gate,speed_gate", [
    ("float", {"float_heave": 5}, 0, .016, .030),
    ("spar", {"spar_heave": -5}, 1, .0025, .0007),
])
def test_published_pto_extension_fine_step_pair(
    name, coordinate, active, position_gate, speed_gate,
):
    hydro = (Path(APPLICATIONS)
             / "_Common_Input_Files/RM3/hydroData/rm3.h5")
    with h5py.File(hydro) as h5:
        center_gap = (np.asarray(h5["body1/properties/cg"]).ravel()[2]
                      - np.asarray(h5["body2/properties/cg"]).ravel()[2])
    source = {}
    python = {}
    for dt, suffix in ((.1, "0p1"), (.05, "0p05"), (.025, "0p025")):
        source[dt] = np.loadtxt(
            Path(SOURCE) / f"{name}_dt_{suffix}.csv", delimiter=",",
        )
        assert source[dt].shape == (round(30 / dt) + 1, 9)
        assert np.isfinite(source[dt]).all()
        python[dt] = _python_response(hydro, coordinate, dt)
        np.testing.assert_allclose(python[dt].time, source[dt][:, 0],
                                   rtol=0, atol=1e-10)

    source_coarse_change = np.max(np.abs(
        source[.1][:, 1:7] - source[.05][::2, 1:7]), axis=0,
    )
    source_fine_change = np.max(np.abs(
        source[.05][:, 1:7] - source[.025][::2, 1:7]), axis=0,
    )
    active_columns = (1, 2, 5, 6) if active == 0 else (3, 4, 5, 6)
    for column in active_columns:
        assert source_fine_change[column - 1] < source_coarse_change[column - 1]

    position_errors = {}
    speed_errors = {}
    for dt in (.05, .025):
        response = python[dt]
        expected = source[dt]
        position_errors[dt] = np.max(np.abs(
            response.body_position[:, active, 2]
            - expected[:, 1 + 2 * active],
        ))
        speed_errors[dt] = np.max(np.abs(
            response.body_velocity[:, active, 2]
            - expected[:, 2 + 2 * active],
        ))
        inactive = 1 - active
        assert np.max(np.abs(response.body_position[:, inactive, 2]
                             - expected[:, 1 + 2 * inactive])) < 1e-8
        assert np.max(np.abs(response.body_velocity[:, inactive, 2]
                             - expected[:, 2 + 2 * inactive])) < 1e-8

        pitch = response.body_position[:, 0, 4]
        pitch_speed = response.body_velocity[:, 0, 4]
        cosine = np.cos(pitch)
        stroke = ((response.body_position[:, 0, 2]
                   - response.body_position[:, 1, 2]) / cosine - center_gap)
        stroke_speed = ((response.body_velocity[:, 0, 2]
                         - response.body_velocity[:, 1, 2]
                         + (center_gap + stroke) * np.sin(pitch) * pitch_speed)
                        / cosine)
        if dt == .025:
            assert np.max(np.abs(stroke - expected[:, 5])) < position_gate
            assert np.max(np.abs(stroke_speed - expected[:, 6])) < speed_gate
        assert np.max(np.abs(response.pto_force - expected[:, 7])) < 1e-7
        assert np.max(np.abs(response.pto_force * stroke_speed
                             - expected[:, 8])) < 1e-7

    assert position_errors[.025] < position_gate
    assert speed_errors[.025] < speed_gate
    assert position_errors[.025] < .75 * position_errors[.05]
    assert speed_errors[.025] < .75 * speed_errors[.05]

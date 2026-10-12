"""Probe the WaveStar plant before NMPC using a derived fine MATLAB step."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim.environment.irregularWave import (
    jonswap_equal_energy_components, synthesize_irregular_response,
)
from wecsim.dynamics.wavestar import run_wavestar_published


REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
HYDRO = os.environ.get("WEC_SIM_WAVESTAR_H5")
pytestmark = pytest.mark.skipif(
    not (REFERENCE and HYDRO),
    reason="derived WaveStar MATLAB output and hydrodynamics are absent",
)


def _source(name):
    values = np.loadtxt(
        Path(REFERENCE) / f"WECCCOMP_NMPC_FINE_DIAG_{name}.csv",
        delimiter=",", ndmin=2,
    )
    assert np.isfinite(values).all()
    return values


def test_fine_step_precontrol_plant_diagnostic():
    components = _source("components")
    wave = _source("wave")
    source_float = _source("WECCCOMP_Nonlinear_Model_Predictive_body1")
    source_pto = _source("WECCCOMP_Nonlinear_Model_Predictive_pto1")
    source_forces = _source("body1_forces")
    assert components.shape == (500, 4)
    assert wave.shape == (9951, 2)
    assert source_float.shape == (9951, 25)
    assert source_pto.shape == (9951, 43)
    assert source_forces.shape == (9951, 37)
    sea = jonswap_equal_energy_components(
        HYDRO, significant_height=.1042, peak_period=1.836,
        directions=np.array([0.]), spreading=np.array([1.]),
        gamma=3.3, phase=components[:, 3, None],
    )
    response = run_wavestar_published(
        HYDRO, sea, dt=.001, end_time=9.95, ramp_time=25,
        g=9.80665,
    )
    incident = synthesize_irregular_response(
        HYDRO, sea, dt=.001, end_time=9.95,
        ramp_time=25, g=9.80665,
    )
    np.testing.assert_allclose(source_float[:, 0], response.time,
                               rtol=0, atol=1e-9)
    np.testing.assert_allclose(wave[:, 1], incident.elevation,
                               rtol=0, atol=1e-10)
    np.testing.assert_allclose(source_float[:, 19:25],
                               response.excitation_force,
                               rtol=0, atol=1e-7)
    assert np.max(np.abs(source_pto[:, 31:37])) < 1e-12
    errors = {
        "pitch_rad": np.max(np.abs(response.angle - source_float[:, 5])),
        "pitch_speed_rad_s": np.max(np.abs(
            response.angular_speed - source_float[:, 11])),
        "float_center_m": np.max(np.abs(
            response.float_position[:, [0, 2]]
            - source_float[:, [1, 3]])),
        "float_center_speed_m_s": np.max(np.abs(
            response.float_velocity[:, [0, 2]]
            - source_float[:, [7, 9]])),
        "pto_stroke_m": np.max(np.abs(
            response.pto_stroke - source_pto[:, 3])),
        "pto_speed_m_s": np.max(np.abs(
            response.pto_speed - source_pto[:, 9])),
        "radiation_N_or_Nm": np.max(np.abs(
            response.radiation_force - source_forces[:, 1:7])),
    }
    limits = {
        "pitch_rad": 2e-5,
        "pitch_speed_rad_s": 4e-4,
        "float_center_m": 1e-5,
        "float_center_speed_m_s": 2e-4,
        "pto_stroke_m": 4e-6,
        "pto_speed_m_s": 8e-5,
        "radiation_N_or_Nm": 5e-4,
    }
    for name, limit in limits.items():
        assert errors[name] < limit, f"{name}: {errors[name]} >= {limit}"
    print("fine-step WaveStar precontrol errors:", errors)

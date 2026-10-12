"""Independent WaveStar fault motion against the pinned MATLAB application."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim.environment.irregularWave import jonswap_equal_energy_components
from wecsim.reference_models.wavestarFault import run_wavestar_fault_published


REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
HYDRO = os.environ.get("WEC_SIM_WAVESTAR_H5")
pytestmark = pytest.mark.skipif(
    not (REFERENCE and HYDRO),
    reason="pinned WaveStar fault output and hydrodynamics are absent",
)


def _source(name):
    values = np.loadtxt(
        Path(REFERENCE) / f"WECCCOMP_FAULT_SOURCE_{name}.csv",
        delimiter=",", ndmin=2,
    )
    assert np.isfinite(values).all()
    return values


def test_independent_published_fault_trajectory():
    components = _source("components")
    true_position = _source("fault_true_position")
    measured_position = _source("fault_measured_position")
    assert true_position.shape == measured_position.shape == (141_201, 2)
    sea = jonswap_equal_energy_components(
        HYDRO, significant_height=.0625, peak_period=1.412,
        directions=np.array([0.]), spreading=np.array([1.]),
        phase=components[:, 3, None],
    )
    response = run_wavestar_fault_published(
        HYDRO, sea,
        measured_position[:, 1] - true_position[:, 1],
        measured_position[:, 1] == 0,
    )
    source_float = _source("WECCCOMP_Fault_Implementation_body1")
    source_arm = _source("WECCCOMP_Fault_Implementation_body2")
    source_pto = _source("WECCCOMP_Fault_Implementation_pto1")
    source_forces = _source("body1_forces")
    source_wave = _source("wave")
    assert response.time.shape == (14_121,)
    np.testing.assert_allclose(response.time, source_float[:, 0],
                               rtol=0, atol=1e-8)
    np.testing.assert_allclose(response.wave_elevation, source_wave[:, 1],
                               rtol=0, atol=1e-10)
    np.testing.assert_allclose(response.excitation_force,
                               source_float[:, 19:25], rtol=0, atol=1e-7)
    np.testing.assert_allclose(response.angle, source_float[:, 5],
                               rtol=0, atol=2e-4)
    np.testing.assert_allclose(response.angular_speed, source_float[:, 11],
                               rtol=0, atol=.015)
    for position, velocity, source in (
        (response.float_position, response.float_velocity, source_float),
        (response.arm_position, response.arm_velocity, source_arm),
    ):
        np.testing.assert_allclose(position[:, [0, 2]], source[:, [1, 3]],
                                   rtol=0, atol=1e-4)
        np.testing.assert_allclose(velocity[:, [0, 2]], source[:, [7, 9]],
                                   rtol=0, atol=.007)
        np.testing.assert_allclose(position[:, 4], source[:, 5],
                                   rtol=0, atol=2e-4)
        np.testing.assert_allclose(velocity[:, 4], source[:, 11],
                                   rtol=0, atol=.015)
    np.testing.assert_allclose(response.pto_stroke, source_pto[:, 3],
                               rtol=0, atol=4e-5)
    np.testing.assert_allclose(response.pto_speed, source_pto[:, 9],
                               rtol=0, atol=.003)
    np.testing.assert_allclose(response.pto_force, source_pto[:, 33],
                               rtol=0, atol=.3)
    np.testing.assert_allclose(response.radiation_force,
                               source_forces[:, 1:7], rtol=0, atol=.01)

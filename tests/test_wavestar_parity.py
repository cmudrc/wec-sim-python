"""Independent WaveStar motion against the pinned published MATLAB case."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim.environment.irregularWave import jonswap_equal_energy_components
from wecsim.reference_models.wavestar import run_wavestar_published


REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
HYDRO = os.environ.get("WEC_SIM_WAVESTAR_H5")
pytestmark = pytest.mark.skipif(
    not (REFERENCE and HYDRO),
    reason="pinned WaveStar output and hydrodynamics are absent",
)


def _source(name):
    values = np.loadtxt(
        Path(REFERENCE) / f"WECCCOMP_SOURCE_{name}.csv",
        delimiter=",", ndmin=2,
    )
    assert np.isfinite(values).all()
    return values


def test_independent_published_wave_star_trajectory():
    source_components = _source("components")
    components = jonswap_equal_energy_components(
        HYDRO, significant_height=.0625, peak_period=1.412,
        directions=np.array([0.]), spreading=np.array([1.]),
        gamma=1, phase=source_components[:, 3, None],
    )
    response = run_wavestar_published(HYDRO, components)
    source_float = _source("WECCCOMP_body1")
    source_arm = _source("WECCCOMP_body2")
    source_pto = _source("WECCCOMP_pto1")
    source_wave = _source("wave")
    source_forces = _source("body1_forces")
    assert response.time.shape == (14_121,)
    np.testing.assert_allclose(response.time, source_float[:, 0],
                               rtol=0, atol=1e-8)
    np.testing.assert_allclose(response.wave_elevation, source_wave[:, 1],
                               rtol=0, atol=1e-10)
    np.testing.assert_allclose(response.excitation_force,
                               source_float[:, 19:25], rtol=0, atol=1e-7)
    for position, velocity, source in (
        (response.float_position, response.float_velocity, source_float),
        (response.arm_position, response.arm_velocity, source_arm),
    ):
        np.testing.assert_allclose(position[:, [0, 2]], source[:, [1, 3]],
                                   rtol=0, atol=1e-5)
        np.testing.assert_allclose(position[:, 4], source[:, 5],
                                   rtol=0, atol=3e-5)
        np.testing.assert_allclose(velocity[:, [0, 2]], source[:, [7, 9]],
                                   rtol=0, atol=1e-4)
        np.testing.assert_allclose(velocity[:, 4], source[:, 11],
                                   rtol=0, atol=2e-4)
    np.testing.assert_allclose(response.pto_stroke, source_pto[:, 3],
                               rtol=0, atol=5e-6)
    np.testing.assert_allclose(response.pto_speed, source_pto[:, 9],
                               rtol=0, atol=5e-5)
    np.testing.assert_allclose(response.radiation_force,
                               source_forces[:, 1:7], rtol=0, atol=1e-3)

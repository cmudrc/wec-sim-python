"""Published fixed hydro monopile against pinned MATLAB WEC-Sim."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim import run_fixed_hydro_monopile
from wecsim.environment.irregularWave import (
    pm_equal_energy_components,
)


REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
HYDRO = os.environ.get("WEC_SIM_MONOPILE_H5")
pytestmark = pytest.mark.skipif(
    not REFERENCE or not HYDRO,
    reason="pinned hydro monopile MATLAB output and HDF5 absent",
)


def _source(name):
    return np.loadtxt(Path(REFERENCE) / f"MONOPILE_HYDRO_{name}.csv",
                      delimiter=",", ndmin=2)


def test_published_fixed_monopile_wave_and_excitation():
    components, wave = _source("components"), _source("wave")
    body1, body2 = _source("monopile_body1"), _source("monopile_body2")
    force1, force2 = _source("body1_forces"), _source("body2_forces")
    joint1, joint2 = _source("constraint1"), _source("constraint2")
    assert components.shape[1] == 4 and len(components) >= 100
    assert wave.shape == (40_001, 2)
    assert body1.shape == body2.shape == (40_001, 25)
    assert force1.shape == force2.shape == (40_001, 37)
    assert joint1.shape == joint2.shape == (40_001, 25)
    assert np.isfinite(components).all()
    assert np.isfinite(wave).all()
    for body in (body1, body2):
        np.testing.assert_allclose(body[:, 0], wave[:, 0],
                                   rtol=0, atol=1e-10)
        assert np.max(np.abs(body[:, 1:7] - body[0, 1:7])) < 1e-6
        assert np.max(np.abs(body[:, 7:13])) < 1e-6

    sea = pm_equal_energy_components(
        HYDRO, significant_height=2, peak_period=5,
        directions=(0,), spreading=(1,), count=len(components),
        phase=components[:, 3:4],
    )
    np.testing.assert_allclose(sea.omega, components[:, 0],
                               rtol=0, atol=1e-12)
    np.testing.assert_allclose(sea.spectral_amplitude, components[:, 1],
                               rtol=0, atol=1e-12)
    np.testing.assert_allclose(sea.d_omega, components[:, 2],
                               rtol=0, atol=1e-12)
    result = run_fixed_hydro_monopile(
        HYDRO, sea, dt=.01, end_time=400, ramp_time=100,
        rho=1025, g=9.81,
    )
    np.testing.assert_allclose(result.time, wave[:, 0],
                               rtol=0, atol=1e-10)
    np.testing.assert_allclose(result.wave_elevation, wave[:, 1],
                               rtol=0, atol=1e-10)
    np.testing.assert_allclose(result.excitation_force,
                               body1[:, 19:25], rtol=0, atol=1e-5)
    for index, (body, force, joint) in enumerate(
        ((body1, force1, joint1), (body2, force2, joint2))
    ):
        np.testing.assert_allclose(result.body_position[:, index],
                                   body[:, 1:7], rtol=0, atol=1e-8)
        np.testing.assert_allclose(result.body_velocity[:, index],
                                   body[:, 7:13], rtol=0, atol=1e-8)
        np.testing.assert_allclose(result.body_force_total[:, index],
                                   body[:, 13:19], rtol=0, atol=1e-5)
        np.testing.assert_allclose(result.joint_force[:, index],
                                   joint[:, 19:25], rtol=0, atol=1e-5)
        # Radiation, added mass, Morison drag, damping, and acceleration
        # vanish for this published stationary case.
        for start in (1, 7, 19, 25, 31):
            assert np.max(np.abs(force[:, start:start + 6])) < 1e-7

    np.testing.assert_allclose(body1[:, 13:19], body1[:, 19:25],
                               rtol=0, atol=1e-7)
    assert np.max(np.abs(body2[:, 19:25])) < 1e-7

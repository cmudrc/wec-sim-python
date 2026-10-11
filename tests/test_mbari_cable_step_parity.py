"""Paired Cable dynamics against a derived finer-step MATLAB source run."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim import run_mbari_cable


REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
HYDRO = os.environ.get("WEC_SIM_CABLE_H5")
pytestmark = pytest.mark.skipif(
    not REFERENCE or not HYDRO,
    reason="pinned MBARI Cable source output and hydrodynamics absent",
)


def _source(name):
    return np.loadtxt(Path(REFERENCE) / f"{name}.csv", delimiter=",")


def test_finer_source_step_and_independent_cable_motion():
    coarse = np.stack([_source(f"CABLE_SOURCE_Cable_body{i}")
                       for i in (1, 2, 3)], axis=1)
    fine = np.stack([_source(f"CABLE_STEP_body{i}")
                     for i in (1, 2, 3)], axis=1)
    finer = np.stack([_source(f"CABLE_STEP_FINER_body{i}")
                      for i in (1, 2, 3)], axis=1)
    coarse_cable = _source("CABLE_SOURCE_cable1")
    fine_cable = _source("CABLE_STEP_cable_force")
    finer_cable = _source("CABLE_STEP_FINER_cable_force")
    coarse_wave = _source("CABLE_SOURCE_wave")
    fine_wave = _source("CABLE_STEP_wave")
    finer_wave = _source("CABLE_STEP_FINER_wave")
    assert coarse.shape == (12_001, 3, 25)
    assert fine.shape == (24_001, 3, 13)
    assert finer.shape == (48_001, 3, 13)
    assert fine_cable.shape == fine_wave.shape == (24_001, 2)
    assert finer_cable.shape == finer_wave.shape == (48_001, 2)
    np.testing.assert_allclose(fine_wave[::2], coarse_wave,
                               rtol=0, atol=1e-10)
    np.testing.assert_allclose(finer_wave[::2], fine_wave,
                               rtol=0, atol=1e-10)
    np.testing.assert_allclose(fine[:, :, 0][::2], coarse[:, :, 0],
                               rtol=0, atol=1e-10)
    np.testing.assert_allclose(finer[:, :, 0][::2], fine[:, :, 0],
                               rtol=0, atol=1e-10)

    # The published source step changes the snap response materially.
    assert np.max(np.abs(fine[::2, :, 1:7] - coarse[:, :, 1:7])) > .03
    assert np.max(np.abs(fine_cable[::2, 1] - coarse_cable[:, 21])) > 30_000
    # The next source refinement is much smaller on the same wave and grid.
    assert np.max(np.abs(finer[::2, :, 1:7] - fine[:, :, 1:7])) < .002
    assert np.max(np.abs(finer_cable[::2, 1] - fine_cable[:, 1])) < 2_500

    result = run_mbari_cable(HYDRO, dt=.005)
    np.testing.assert_allclose(result.time, fine_wave[:, 0],
                               rtol=0, atol=1e-10)
    np.testing.assert_allclose(result.wave_elevation, fine_wave[:, 1],
                               rtol=0, atol=1e-11)
    position_error = result.body_position - fine[:, :, 1:7]
    speed_error = result.body_velocity - fine[:, :, 7:13]
    force_error = result.cable_force_z - fine_cable[:, 1]
    assert np.max(np.abs(position_error)) < .002
    assert np.max(np.abs(speed_error)) < .025
    assert np.max(np.abs(force_error)) < 2_500
    assert np.sqrt(np.mean(force_error**2)) < 400
    assert abs(np.trapezoid(force_error, result.time)) < 100

    finer_position_error = result.body_position - finer[::2, :, 1:7]
    finer_speed_error = result.body_velocity - finer[::2, :, 7:13]
    finer_force_error = result.cable_force_z - finer_cable[::2, 1]
    assert np.max(np.abs(finer_position_error)) < .001
    assert np.max(np.abs(finer_speed_error)) < .01
    assert np.max(np.abs(finer_force_error)) < 1_000
    assert np.sqrt(np.mean(finer_force_error**2)) < 150
    assert abs(np.trapezoid(finer_force_error, result.time)) < 50

    coarse_python = run_mbari_cable(HYDRO, dt=.01)
    assert np.max(np.abs(coarse_python.body_position
                         - result.body_position[::2])) < 2e-5
    assert np.max(np.abs(coarse_python.body_velocity
                         - result.body_velocity[::2])) < 2.5e-4
    assert np.max(np.abs(coarse_python.cable_force_z
                         - result.cable_force_z[::2])) < 25

"""Check the suspect fitted radiation against the supplied source data."""

import os
from pathlib import Path

import numpy as np
import pytest

from tools.rm3_radiation_diagnostic import (
    active_mode_minimum_damping, common_surge_curves,
    common_surge_feedthrough,
)
from wecsim.hydrodynamics.radiation import RadiationStateSpace
from wecsim.dynamics.rm3Regular import solve_rm3_regular


RM3 = (Path(__file__).parent / "test_objects" / "test_bodyclass"
       / "testData" / "hydroData" / "rm3.h5")
CURRENT_RM3 = os.environ.get("WEC_SIM_RM3_H5")


@pytest.mark.skipif(not CURRENT_RM3, reason="current BEMIO RM3 HDF5 unavailable")
def test_pinned_bemio_radiation_transfer_matches_diagnostic_fixture():
    """Check the transfer, since equivalent state-space realizations may differ."""
    for diagnostic in (common_surge_curves, active_mode_minimum_damping):
        fixture = diagnostic(RM3)
        current = diagnostic(CURRENT_RM3)
        for expected, actual in zip(fixture, current):
            np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-6)
    fixture = active_mode_minimum_damping(RM3, return_modes=True)
    current = active_mode_minimum_damping(CURRENT_RM3, return_modes=True)
    wave_index = np.argmin(abs(fixture[2] - 2 * np.pi / 8))
    pitch_index = np.argmin(abs(fixture[2] - 1.2))
    for index in (0, wave_index, pitch_index):
        np.testing.assert_allclose(
            np.abs(current[4][:, index]), np.abs(fixture[4][:, index]),
            rtol=0, atol=1e-5,
        )


def test_fitted_common_surge_has_opposite_low_frequency_sign():
    source_frequency, source, fitted_frequency, fitted = common_surge_curves(RM3)
    feedthrough = common_surge_feedthrough(RM3)
    assert fitted_frequency[0] == 0
    np.testing.assert_array_equal(fitted_frequency[1:], source_frequency)
    assert np.all(source[:, source_frequency < 0.3] >= 0)
    assert np.all(fitted[:, 0] < -10_000)
    # The saved direct terms are omitted by the pinned MATLAB body class.
    # Even if applied, they do not repair the negative low-frequency fit.
    assert np.all(feedthrough > 0)
    assert np.all(fitted[:, 0] + feedthrough < -2_000)
    assert np.all(fitted[:, np.argmin(abs(fitted_frequency - 2 * np.pi / 8))] > 0)


@pytest.mark.parametrize("body_to_body", [False, True])
def test_active_state_space_retains_source_common_surge_transfer(body_to_body):
    """Removing inactive joint channels must preserve the fitted slow mode."""
    model = RadiationStateSpace.from_rm3_hdf5(
        RM3, body_to_body=body_to_body,
    )
    shared_surge_speed = np.zeros(12)
    shared_surge_speed[[0, 6]] = 1
    steady_force = model.C @ np.linalg.solve(
        -model.A, model.B @ shared_surge_speed,
    )
    fitted_dc = common_surge_curves(RM3)[3][int(body_to_body), 0]
    np.testing.assert_allclose(steady_force[[0, 6]].sum(), fitted_dc,
                               rtol=0, atol=1e-8)


def test_fitted_radiation_has_large_negative_active_joint_modes():
    """Compare physical RM3 joint motions, not unconstrained six-DOF bodies."""
    frequency, source, fitted_frequency, fitted = active_mode_minimum_damping(RM3)
    np.testing.assert_array_equal(fitted_frequency[1:], frequency)
    # Source BEM itself has small negative eigenvalues, so this does not
    # certify global BEM passivity. The fitted negative directions are much
    # larger, even in the four-dimensional moving-joint subspace.
    assert np.all(source.min(axis=1) > -500)
    assert fitted[0, 0] < -14_000
    assert fitted[1, 0] < -17_000
    wave_index = np.argmin(abs(frequency - 2 * np.pi / 8))
    assert -50 < source[1, wave_index] < 0
    assert fitted[1, wave_index + 1] < -5_000


def test_fitted_negative_modes_require_a_coupled_refit():
    """Locate three nonpassive directions before considering a physical refit."""
    _, _, fitted_frequency, minimum, modes = active_mode_minimum_damping(
        RM3, return_modes=True,
    )
    for coupled in (False, True):
        dc = np.abs(modes[int(coupled), 0])
        assert minimum[int(coupled), 0] < -14_000
        assert dc[0] > .99  # Common surge dominates near zero frequency.
        assert np.linalg.norm(dc[1:3]) < 1e-3

    wave_index = np.argmin(abs(fitted_frequency - 2 * np.pi / 8))
    wave = np.abs(modes[1, wave_index])
    assert minimum[1, wave_index] < -5_000
    assert wave[0] < 1e-3 and wave[3] < 1e-3
    assert wave[1] > .4 and wave[2] > .9  # Coupled heave dominates at 8 s.

    pitch_index = np.argmin(abs(fitted_frequency - 1.2))
    pitch = np.abs(modes[1, pitch_index])
    assert minimum[1, pitch_index] < -4_000
    assert pitch[0] < .02 and np.linalg.norm(pitch[1:3]) < 1e-3
    assert pitch[3] > .99  # Shared pitch dominates near 1.2 rad/s.


@pytest.mark.parametrize("body_to_body", [False, True])
def test_default_convolution_does_not_amplify_no_wave_common_surge(body_to_body):
    """The physical default should dissipate a free surge perturbation.

    This checks an unforced response independently of the published MATLAB
    Cases 5–6, whose fitted radiation supplies negative low-frequency damping.
    """
    initial_speed = 0.01
    response = solve_rm3_regular(
        RM3, wave_height=0, no_wave=True, b2b=body_to_body,
        radiation_memory=60, end_time=400,
        initial_speed=np.array([initial_speed, 0, 0, 0]),
    )
    surge_speed = response.body_velocity[:, 0, 0]
    np.testing.assert_allclose(surge_speed[0], initial_speed, rtol=0, atol=1e-12)
    assert np.max(np.abs(surge_speed)) <= initial_speed * (1 + 1e-6)
    assert 0 < surge_speed[-1] < initial_speed

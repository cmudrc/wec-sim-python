"""Independently advance the published WaveStar plant before NMPC starts."""

from pathlib import Path

import numpy as np
import pytest

from wecsim import WaveStarNmpcPTO
from wecsim.dynamics.wavestar import run_wavestar_published
from wecsim.environment.irregularWave import jonswap_equal_energy_components


ROOT = Path(__file__).parent / "fixtures/matlab"
HYDRO = ROOT / "wavestar_fixed_input_step/wavestar.h5"
SOURCE = ROOT / "wavestar_fixed_input_step/source.npz"
FULL_SOURCE = ROOT / "wavestar_nmpc_full_step/baseline.npz"


def _published_sea(components):
    return jonswap_equal_energy_components(
        HYDRO, significant_height=.1042, peak_period=1.836,
        directions=np.array([0.]), spreading=np.array([1.]),
        gamma=3.3, phase=components[:, 3, None],
    )


def test_public_source_ode8_closes_the_precontrol_motion_gap():
    with np.load(SOURCE, allow_pickle=False) as published:
        body = published["50ms_WECCCOMP_Nonlinear_Model_Predictive_body1"]
        sea = _published_sea(published["50ms_components"])
    response = run_wavestar_published(
        HYDRO, sea, dt=.05, end_time=9.95, ramp_time=25,
        g=9.80665, solver="wecsim_ode8",
    )
    np.testing.assert_allclose(response.time, body[:, 0], rtol=0, atol=1e-10)
    # Source motion is never fed to the Python solve. Only pinned HDF5
    # hydrodynamics and the published wave phases set its inputs.
    assert np.max(np.abs(response.angle - body[:, 5])) < 2e-5
    assert np.max(np.abs(response.angular_speed - body[:, 11])) < 1.5e-4


def test_wavestar_rejects_unknown_solver():
    with np.load(SOURCE, allow_pickle=False) as published:
        sea = _published_sea(published["50ms_components"])
    with pytest.raises(ValueError, match="WaveStar solver"):
        run_wavestar_published(HYDRO, sea, solver="unavailable")


def test_source_ode8_with_independent_nmpc_through_20s():
    with np.load(FULL_SOURCE, allow_pickle=False) as published:
        body = published["body"][:401]
        pto = published["pto"][:401]
        sea = _published_sea(published["components"])
    response = run_wavestar_published(
        HYDRO, sea, dt=.05, end_time=20, ramp_time=25,
        g=9.80665, solver="wecsim_ode8",
        pto_controller=WaveStarNmpcPTO(plant_dt=.05, control_dt=.05),
    )
    np.testing.assert_allclose(response.time, body[:, 0], rtol=0, atol=1e-10)
    # NMPC activates at 15 s. This bounded diagnostic does not gate the
    # chaotic full 225 s closed-loop motion.
    assert np.max(np.abs(response.angle - body[:, 5])) < .002
    assert np.max(np.abs(response.angular_speed - body[:, 11])) < .012
    assert np.max(np.abs(response.pto_force - pto[:, 33])) < 1.

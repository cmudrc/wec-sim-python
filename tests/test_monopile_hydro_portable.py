"""Full published monopile hydro case, using pinned MATLAB data in Python CI."""

from pathlib import Path

import numpy as np

from wecsim import run_fixed_hydro_monopile
from wecsim.environment.irregularWave import pm_equal_energy_components


FIXTURES = Path(__file__).parent / "fixtures/matlab"


def test_published_fixed_monopile_full_duration():
    """Advance the 400 s sea from phases, then check both fixed-body loads."""
    with np.load(FIXTURES / "monopile_hydro.npz", allow_pickle=False) as source:
        components = source["components"]
        wave = source["wave"]
        assert components.shape == (500, 4)
        assert wave.shape == (40_001, 2)
        assert source["body_position"].shape == (40_001, 2, 6)
        assert source["body_velocity"].shape == (40_001, 2, 6)
        assert source["body_force_total"].shape == (40_001, 2, 6)
        assert source["joint_force"].shape == (40_001, 2, 6)

        hydro = FIXTURES / "monopile_hydro.h5"
        sea = pm_equal_energy_components(
            hydro, significant_height=2, peak_period=5,
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
            hydro, sea, dt=.01, end_time=400, ramp_time=100,
            rho=1025, g=9.81,
        )
        np.testing.assert_allclose(result.time, wave[:, 0], rtol=0,
                                   atol=1e-10)
        np.testing.assert_allclose(result.wave_elevation, wave[:, 1],
                                   rtol=0, atol=1e-10)
        np.testing.assert_allclose(result.excitation_force,
                                   source["body_excitation"], rtol=0,
                                   atol=1e-5)
        np.testing.assert_allclose(result.body_position,
                                   source["body_position"], rtol=0,
                                   atol=1e-8)
        np.testing.assert_allclose(result.body_velocity,
                                   source["body_velocity"], rtol=0,
                                   atol=1e-8)
        np.testing.assert_allclose(result.body_force_total,
                                   source["body_force_total"], rtol=0,
                                   atol=1e-5)
        np.testing.assert_allclose(result.joint_force,
                                   source["joint_force"], rtol=0,
                                   atol=1e-5)

        # The second body's only load is its weight; its joint cancels it.
        tower = result.body_force_total[:, 1]
        np.testing.assert_allclose(tower[:, 2], -1_031_930 * 9.81,
                                   rtol=0, atol=1e-8)
        np.testing.assert_allclose(result.joint_force[:, 1] + tower, 0,
                                   rtol=0, atol=1e-8)

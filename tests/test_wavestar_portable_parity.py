"""Pair the published WaveStar base motion in ordinary Python CI."""

from pathlib import Path

import numpy as np

from wecsim.dynamics.wavestar import run_wavestar_published
from wecsim.environment.irregularWave import jonswap_equal_energy_components


ROOT = Path(__file__).resolve().parents[1]
HYDRO = ROOT / "examples/data/wavestar.h5"
MATLAB = ROOT / "tests/fixtures/matlab/wavestar_base.npz"


def test_full_duration_published_wavestar_base_case():
    """Advance motion from pinned HDF5 and phases, never source state or force."""
    with np.load(MATLAB, allow_pickle=False) as source:
        assert source["time"].shape == (14_121,)
        assert source["components"].shape == (500, 4)
        components = jonswap_equal_energy_components(
            HYDRO, significant_height=.0625, peak_period=1.412,
            directions=np.array([0.]), spreading=np.array([1.]),
            gamma=1, phase=source["components"][:, 3, None],
        )
        for actual, column in (
            (components.omega, 0),
            (components.spectral_amplitude, 1),
            (components.d_omega, 2),
        ):
            np.testing.assert_allclose(actual, source["components"][:, column],
                                       rtol=0, atol=1e-12)

        result = run_wavestar_published(HYDRO, components)
        np.testing.assert_allclose(result.time, source["time"],
                                   rtol=0, atol=1e-8)
        np.testing.assert_allclose(result.wave_elevation,
                                   source["wave_elevation"],
                                   rtol=0, atol=1e-10)
        np.testing.assert_allclose(result.excitation_force,
                                   source["float_excitation"],
                                   rtol=0, atol=1e-7)
        for position, velocity, label in (
            (result.float_position, result.float_velocity, "float"),
            (result.arm_position, result.arm_velocity, "arm"),
        ):
            np.testing.assert_allclose(position[:, [0, 2]],
                                       source[f"{label}_position"][:, [0, 2]],
                                       rtol=0, atol=1e-5)
            np.testing.assert_allclose(position[:, 4],
                                       source[f"{label}_position"][:, 4],
                                       rtol=0, atol=3e-5)
            np.testing.assert_allclose(velocity[:, [0, 2]],
                                       source[f"{label}_velocity"][:, [0, 2]],
                                       rtol=0, atol=1e-4)
            np.testing.assert_allclose(velocity[:, 4],
                                       source[f"{label}_velocity"][:, 4],
                                       rtol=0, atol=2e-4)
            for inactive in (1, 3, 5):
                np.testing.assert_allclose(position[:, inactive],
                                           source[f"{label}_position"][:, inactive],
                                           rtol=0, atol=1e-10)
                np.testing.assert_allclose(velocity[:, inactive],
                                           source[f"{label}_velocity"][:, inactive],
                                           rtol=0, atol=1e-10)

        np.testing.assert_allclose(result.pto_stroke, source["pto_stroke"],
                                   rtol=0, atol=5e-6)
        np.testing.assert_allclose(result.pto_speed, source["pto_speed"],
                                   rtol=0, atol=5e-5)
        np.testing.assert_allclose(result.radiation_force,
                                   source["radiation_force"],
                                   rtol=0, atol=1e-3)
        np.testing.assert_array_equal(result.pto_force, 0)
        assert np.max(np.abs(result.angle)) > .1

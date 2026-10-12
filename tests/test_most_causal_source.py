"""Pair causal MOST coupling with the full published source trajectory."""

import os
from pathlib import Path

import numpy as np
import pytest
from scipy.io import loadmat

from wecsim import (
    MostBaselineController, MostCoupled, MostPlatformHydrodynamics,
    MostRotor, MostTowerReaction, MostWindField, read_turbsim_bts,
)
from wecsim.environment.irregularWave import (
    jonswap_equal_energy_components, synthesize_irregular_response,
)


@pytest.mark.skipif(
    not all(os.environ.get(name) for name in (
        "WEC_SIM_MOST_FULL_SOURCE", "WEC_SIM_MOST_H5",
        "WEC_SIM_MOST_MASS_PROPERTIES", "WEC_SIM_MOST_PROPERTIES",
        "WEC_SIM_MOST_CONTROL", "WEC_SIM_MOST_STEADY_STATES",
        "WEC_SIM_MOST_BLADE_DIR", "WEC_SIM_MOST_BTS",
    )),
    reason="pinned published MOST trajectory and generated inputs not provided",
)
def test_most_causal_1000s_against_published_source():
    hydro = os.environ["WEC_SIM_MOST_H5"]
    platform = MostPlatformHydrodynamics.from_volturnus(
        hydro, os.environ["WEC_SIM_MOST_MASS_PROPERTIES"],
    )
    tower = MostTowerReaction.from_iea15mw(
        os.environ["WEC_SIM_MOST_PROPERTIES"],
        platform_cg=platform.equilibrium_pose[:3],
    )
    controller = MostBaselineController.from_matlab_files(
        os.environ["WEC_SIM_MOST_CONTROL"],
        os.environ["WEC_SIM_MOST_STEADY_STATES"], wind_speed=8,
    )
    rotor = MostRotor.from_iea15mw(
        Path(os.environ["WEC_SIM_MOST_BLADE_DIR"]), controller=controller,
    )
    wind = MostWindField(read_turbsim_bts(os.environ["WEC_SIM_MOST_BTS"]))
    sea = jonswap_equal_energy_components(
        hydro, significant_height=4, peak_period=8,
        directions=np.array([0.]), spreading=np.array([1.]),
        seed=1, phase_generator="matlab",
    )
    wave = synthesize_irregular_response(
        hydro, sea, dt=.01, end_time=1000, ramp_time=20,
        rho=1025, g=9.80665,
    )
    assert wave.time.size == 100001
    result = MostCoupled(platform, rotor, tower).simulate_causal(
        wave.time, wave.excitation_force, wind,
    )
    assert result.iterations == 1
    source = loadmat(os.environ["WEC_SIM_MOST_FULL_SOURCE"])
    np.testing.assert_allclose(
        result.platform.time[::10], source["body_time"][:, 0],
        rtol=0, atol=1e-12,
    )
    for axis, position_gate, velocity_gate in (
        (0, 2e-3, 1.5e-4),
        (1, 11e-3, 7e-4),
        (2, 2e-4, 5e-5),
        (3, 2.3e-4, 5e-5),
        (4, 5.5e-5, 1.2e-5),
        (5, 4e-4, 4e-5),
    ):
        np.testing.assert_allclose(
            result.platform.position[::10, axis],
            source["body_position"][:, axis],
            rtol=0, atol=position_gate,
        )
        np.testing.assert_allclose(
            result.platform.velocity[::10, axis],
            source["body_velocity"][:, axis],
            rtol=0, atol=velocity_gate,
        )
    np.testing.assert_allclose(
        result.rotor.rotor_speed[::10] * 60 / (2 * np.pi),
        source["rotor_speed"][:, 0], rtol=0, atol=1.2e-3,
    )
    np.testing.assert_allclose(
        result.rotor.azimuth[::10], source["azimuth"][:, 0],
        rtol=0, atol=2e-3,
    )
    np.testing.assert_allclose(
        result.rotor.generator_torque[::10],
        source["generator_torque"][:, 0], rtol=0, atol=5000,
    )
    expected_load = source["blade_aero_load"]
    component_peak = np.max(np.abs(expected_load), axis=(0, 2))
    component_error = np.max(
        np.abs(result.rotor.blade_root_load[::10] - expected_load),
        axis=(0, 2),
    )
    assert np.all(component_error < .04 * component_peak), (
        component_error, component_peak,
    )

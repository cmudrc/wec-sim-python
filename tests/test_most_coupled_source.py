"""Pair a source-free six-DOF MOST platform/turbine trajectory with MATLAB."""

import os
from pathlib import Path

import numpy as np
import pytest
from scipy.io import loadmat

from wecsim import (
    MostBaselineController, MostCoupled, MostPlatformHydrodynamics, MostRotor,
    MostStaticMooring, MostTowerReaction, MostWindField, read_turbsim_bts,
)
from wecsim.environment.irregularWave import (
    jonswap_equal_energy_components, synthesize_irregular_response,
)


@pytest.mark.skipif(
    not all(os.environ.get(name) for name in (
        "WEC_SIM_MOST_SHORT_BASELINE", "WEC_SIM_MOST_H5",
        "WEC_SIM_MOST_MASS_PROPERTIES", "WEC_SIM_MOST_PROPERTIES",
        "WEC_SIM_MOST_BLADE_DIR", "WEC_SIM_MOST_ADVECTION_DIR",
    )),
    reason="pinned MATLAB MOST coupled-run inputs not provided",
)
def test_most_six_dof_coupled_trajectory_against_pinned_source():
    h5_file = os.environ["WEC_SIM_MOST_H5"]
    platform = MostPlatformHydrodynamics.from_volturnus(
        h5_file, os.environ["WEC_SIM_MOST_MASS_PROPERTIES"],
    )
    tower = MostTowerReaction.from_iea15mw(
        os.environ["WEC_SIM_MOST_PROPERTIES"],
        platform_cg=platform.equilibrium_pose[:3],
    )
    controller = None
    if os.environ.get("WEC_SIM_MOST_CONTROL"):
        controller = MostBaselineController.from_matlab_files(
            os.environ["WEC_SIM_MOST_CONTROL"],
            os.environ["WEC_SIM_MOST_STEADY_STATES"], wind_speed=8,
        )
    rotor = MostRotor.from_iea15mw(
        os.environ["WEC_SIM_MOST_BLADE_DIR"], controller=controller,
    )
    wind = MostWindField(read_turbsim_bts(
        Path(os.environ["WEC_SIM_MOST_ADVECTION_DIR"]) / "WIND_8mps.bts",
    ))
    sea = jonswap_equal_energy_components(
        h5_file,
        significant_height=float(os.environ.get("WEC_SIM_MOST_WAVE_HEIGHT", "4")),
        peak_period=8,
        directions=np.array([0.]), spreading=np.array([1.]),
        seed=int(os.environ.get("WEC_SIM_MOST_PHASE_SEED", "1")),
        phase_generator="matlab",
    )
    end_time = float(os.environ.get("WEC_SIM_MOST_END_TIME", "10"))
    assert end_time in (10, 30, 60)
    developed_sea = end_time >= 30
    extended_sea = end_time == 60
    if developed_sea:
        assert controller is not None, "developed-sea case needs its generated MATLAB controller"
    wave = synthesize_irregular_response(
        h5_file, sea, dt=.01, end_time=end_time, ramp_time=20,
        rho=1025, g=9.80665,
    )
    time = wave.time
    assert time.shape == (round(end_time/.01) + 1,)

    # The source trajectory is used only after the independent solve.
    result = MostCoupled(platform, rotor, tower).simulate(
        time, wave.excitation_force, wind,
    )
    source = loadmat(os.environ["WEC_SIM_MOST_SHORT_BASELINE"])
    if controller is not None:
        np.testing.assert_allclose(
            controller.initial_omega*60/(2*np.pi), source["rotor_speed"][0, 0],
            rtol=0, atol=1e-12,
        )
        initial_state = np.array([0, controller.initial_omega/5, 0])
        np.testing.assert_allclose(
            controller._command(initial_state, controller.initial_pitch)[0],
            source["generator_torque"][0, 0], rtol=0, atol=1e-6,
        )
    np.testing.assert_allclose(time, source["body_time"].ravel(),
                               rtol=0, atol=1e-12)
    for actual, expected in (
        (sea.omega, source["wave_omega"].ravel()),
        (sea.spectral_amplitude, source["wave_amplitude"].ravel()),
        (sea.d_omega, source["wave_d_omega"].ravel()),
        (sea.phase, source["wave_phase"]),
    ):
        np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-11)
    np.testing.assert_allclose(
        wave.elevation, source["wave_elevation"].ravel(), rtol=0, atol=1e-12,
    )
    np.testing.assert_allclose(
        wave.excitation_force, source["body_force_excitation"],
        rtol=0, atol=1e-4,
    )
    assert result.iterations <= (20 if developed_sea else 12)
    assert result.position_residual <= 1e-6
    assert result.velocity_residual <= 1e-6
    for axis, position_gate, velocity_gate in (
        (0, 6e-4 if developed_sea else 2e-4, 5e-5),  # surge, m and m/s
        (1, 7.5e-4 if extended_sea else 5e-4, 2e-4),  # sway, m and m/s
        (2, 1.5e-4, 7e-5),  # heave, m and m/s
        (3, 6e-5 if extended_sea else 4e-5, 2e-5),  # roll, rad and rad/s
        (4, 8e-6 if extended_sea else 5e-6, 2e-6),  # pitch, rad and rad/s
        (5, 1.5e-4, 1e-5 if developed_sea else 3e-6),  # yaw, rad and rad/s
    ):
        np.testing.assert_allclose(
            result.platform.position[:, axis], source["body_position"][:, axis],
            rtol=0, atol=position_gate,
        )
        np.testing.assert_allclose(
            result.platform.velocity[:, axis], source["body_velocity"][:, axis],
            rtol=0, atol=velocity_gate,
        )
    np.testing.assert_allclose(
        result.rotor.rotor_speed*60/(2*np.pi), source["rotor_speed"].ravel(),
        rtol=0, atol=.001,
    )
    np.testing.assert_allclose(
        result.rotor.azimuth, source["azimuth"].ravel(),
        rtol=0, atol=.0004 if developed_sea else .0003,
    )
    np.testing.assert_allclose(
        result.rotor.generator_torque, source["generator_torque"].ravel(),
        rtol=0, atol=1000,
    )
    if developed_sea:
        expected_load = source["blade_aero_load"]
        assert result.rotor.blade_root_load.shape == expected_load.shape == (
            len(time), 6, 3,
        )
        component_peak = np.max(np.abs(expected_load), axis=(0, 2))
        component_error = np.max(
            np.abs(result.rotor.blade_root_load - expected_load), axis=(0, 2),
        )
        assert np.all(component_error < 0.04 * component_peak), (
            component_error, component_peak,
        )


@pytest.mark.skipif(
    not all(os.environ.get(name) for name in (
        "WEC_SIM_MOST_SHORT_BASELINE", "WEC_SIM_MOST_H5",
        "WEC_SIM_MOST_MASS_PROPERTIES", "WEC_SIM_MOST_PROPERTIES",
    )),
    reason="pinned MATLAB MOST force history not provided",
)
def test_most_source_motion_force_laws_against_pinned_source():
    """Replay 30/60 s force laws on saved source motion, independent of solve."""
    source = loadmat(os.environ["WEC_SIM_MOST_SHORT_BASELINE"])
    platform = MostPlatformHydrodynamics.from_volturnus(
        os.environ["WEC_SIM_MOST_H5"],
        os.environ["WEC_SIM_MOST_MASS_PROPERTIES"],
    )
    np.testing.assert_allclose(
        platform.restoring_force(source["body_position"]),
        source["body_force_restoring"], rtol=0, atol=1e-5,
    )
    np.testing.assert_allclose(
        platform.drag_force(source["body_velocity"]),
        source["body_force_viscous"], rtol=0, atol=1e-7,
    )
    np.testing.assert_allclose(
        platform.radiation_force(source["body_velocity"], .01),
        source["body_force_radiation"], rtol=0, atol=1e-6,
    )
    mooring = MostStaticMooring()
    mooring_force = np.array([
        mooring.force(pose)[0] for pose in source["mooring_position"]
    ])
    np.testing.assert_allclose(
        mooring_force[:, :3], source["mooring_force"][:, :3],
        rtol=0, atol=1e-3,
    )
    np.testing.assert_allclose(
        mooring_force[:, 3:], source["mooring_force"][:, 3:],
        rtol=0, atol=1e-2,
    )
    tower = MostTowerReaction.from_iea15mw(
        os.environ["WEC_SIM_MOST_PROPERTIES"],
        platform_cg=platform.equilibrium_pose[:3],
    )
    tower_force = tower.evaluate(
        source["body_position"], source["body_velocity"],
        source["body_acceleration"],
        source["rotor_speed"].ravel()*2*np.pi/60,
        source["azimuth"].ravel(), source["generator_torque"].ravel(),
        source["blade_aero_load"],
    )
    np.testing.assert_allclose(
        tower_force[:, :3], source["tower_base_load"][:, :3],
        rtol=0, atol=1e-5,
    )
    np.testing.assert_allclose(
        tower_force[:, 3:], source["tower_base_load"][:, 3:],
        rtol=0, atol=0.1,
    )

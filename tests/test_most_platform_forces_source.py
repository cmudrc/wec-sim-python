"""Pair active MOST wave and mooring forces on a moving platform."""

import os

import numpy as np
import pytest
from scipy.io import loadmat

from wecsim import MostPlatformHydrodynamics, MostStaticMooring, MostTowerReaction
from wecsim.environment.irregularWave import (
    jonswap_equal_energy_components, synthesize_irregular_response,
)


@pytest.mark.skipif(
    not all(os.environ.get(name) for name in (
        "WEC_SIM_MOST_SHORT_BASELINE", "WEC_SIM_MOST_H5",
        "WEC_SIM_MOST_MASS_PROPERTIES",
    )),
    reason="pinned MATLAB MOST platform force trace not provided",
)
def test_most_platform_wave_and_mooring_against_pinned_source():
    source = loadmat(os.environ["WEC_SIM_MOST_SHORT_BASELINE"])
    h5_file = os.environ["WEC_SIM_MOST_H5"]
    time = source["body_time"].ravel()
    np.testing.assert_allclose(time, source["wave_time"].ravel(),
                               rtol=0, atol=1e-12)
    assert time.shape == (1001,)

    sea = jonswap_equal_energy_components(
        h5_file, significant_height=4, peak_period=8,
        directions=np.array([0.]), spreading=np.array([1.]),
        seed=1, phase_generator="matlab",
    )
    for actual, expected in (
        (sea.omega, source["wave_omega"].ravel()),
        (sea.spectral_amplitude, source["wave_amplitude"].ravel()),
        (sea.d_omega, source["wave_d_omega"].ravel()),
        (sea.phase, source["wave_phase"]),
    ):
        np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-12)
    wave = synthesize_irregular_response(
        h5_file, sea, dt=.01, end_time=10, ramp_time=20,
        rho=1025, g=9.80665,
    )
    np.testing.assert_allclose(wave.elevation,
                               source["wave_elevation"].ravel(),
                               rtol=0, atol=1e-11)
    np.testing.assert_allclose(wave.excitation_force,
                               source["body_force_excitation"],
                               rtol=0, atol=1e-5)

    mooring = MostStaticMooring()
    force = np.array([
        mooring.force(pose)[0] for pose in source["mooring_position"]
    ])
    np.testing.assert_allclose(force[:, :3],
                               source["mooring_force"][:, :3],
                               rtol=0, atol=1e-3)
    np.testing.assert_allclose(force[:, 3:],
                               source["mooring_force"][:, 3:],
                               rtol=0, atol=1e-2)

    platform = MostPlatformHydrodynamics.from_volturnus(
        h5_file, os.environ["WEC_SIM_MOST_MASS_PROPERTIES"],
    )
    np.testing.assert_allclose(
        platform.restoring_force(source["body_position"]),
        source["body_force_restoring"], rtol=0, atol=1e-5,
    )
    np.testing.assert_allclose(
        platform.drag_force(source["body_velocity"]),
        source["body_force_viscous"], rtol=0, atol=1e-8,
    )
    np.testing.assert_allclose(
        platform.radiation_force(source["body_velocity"], .01),
        source["body_force_radiation"], rtol=0, atol=1e-6,
    )

    # Logged WEC-Sim body-force channels use resisting-force signs.
    total = (source["body_force_excitation"]
             - source["body_force_radiation"]
             - source["body_force_added_mass"]
             - source["body_force_restoring"]
             - source["body_force_viscous"]
             - source["body_force_linear_damping"])
    np.testing.assert_allclose(total, source["body_force_total"],
                               rtol=0, atol=1e-5)
    assert np.ptp(source["mooring_position"][:, 0]) > 1

    # The turbine tower-base load is logged in the rotating platform frame.
    # WEC-Sim moves twice the trace of translational infinite-frequency
    # added mass into the Simscape body's scalar mass in this case.
    adjusted_mass = platform.mass + 2*np.trace(platform.added_mass[:3, :3])
    rotation = np.array([
        mooring._rotation(*pose[3:]) for pose in source["body_position"]
    ])
    tower_force_world = np.einsum(
        "nij,nj->ni", rotation, source["tower_base_load"][:, :3],
    )
    net_force = (source["body_force_total"][:, :3]
                 + source["mooring_force"][:, :3] + tower_force_world)
    np.testing.assert_allclose(
        net_force, adjusted_mass*source["body_acceleration"][:, :3],
        rtol=0, atol=1e-4,
    )

    # Advance Python platform motion from Python wave/hydro/mooring. The
    # source tower reaction is the sole prescribed dynamic load; only the
    # published dominant surge, heave, and pitch coordinates are advanced.
    response = platform.simulate(
        time, wave.excitation_force, source["tower_base_load"],
        mooring=mooring,
    )
    for column, position_gate, velocity_gate in (
        (0, 3e-4, 7e-5),  # surge, m and m/s
        (2, 1.5e-4, 7e-5),  # heave, m and m/s
        (4, 1e-5, 1e-5),  # pitch, rad and rad/s
    ):
        np.testing.assert_allclose(
            response.position[:, column], source["body_position"][:, column],
            rtol=0, atol=position_gate,
        )
        np.testing.assert_allclose(
            response.velocity[:, column], source["body_velocity"][:, column],
            rtol=0, atol=velocity_gate,
        )


@pytest.mark.skipif(
    not all(os.environ.get(name) for name in (
        "WEC_SIM_MOST_SHORT_BASELINE", "WEC_SIM_MOST_PROPERTIES",
    )),
    reason="pinned MATLAB MOST turbine and tower trace not provided",
)
def test_most_tower_reaction_against_pinned_source():
    source = loadmat(os.environ["WEC_SIM_MOST_SHORT_BASELINE"])
    assert source["body_position"].shape == (1001, 6)
    tower = MostTowerReaction.from_iea15mw(os.environ["WEC_SIM_MOST_PROPERTIES"])
    force = tower.evaluate(
        source["body_position"], source["body_velocity"],
        source["body_acceleration"],
        source["rotor_speed"].ravel()*2*np.pi/60,
        source["azimuth"].ravel(), source["generator_torque"].ravel(),
        source["blade_aero_load"],
    )
    source_force = source["tower_base_load"]
    assert source_force.shape == (1001, 6)
    np.testing.assert_allclose(force[:, :3], source_force[:, :3],
                               rtol=0, atol=1e-5)
    assert np.all(np.max(np.abs(force[:, 3:] - source_force[:, 3:]), axis=0)
                  < [0.1, 1e-4, 0.01])


@pytest.mark.skipif(
    not all(os.environ.get(name) for name in (
        "WEC_SIM_MOST_SHORT_BASELINE", "WEC_SIM_MOST_H5",
        "WEC_SIM_MOST_MASS_PROPERTIES", "WEC_SIM_MOST_PROPERTIES",
    )),
    reason="pinned MATLAB MOST coupled-platform inputs not provided",
)
def test_most_platform_six_dof_with_implicit_turbine_against_pinned_source():
    source = loadmat(os.environ["WEC_SIM_MOST_SHORT_BASELINE"])
    time = source["body_time"].ravel()
    h5_file = os.environ["WEC_SIM_MOST_H5"]
    sea = jonswap_equal_energy_components(
        h5_file, significant_height=4, peak_period=8,
        directions=np.array([0.]), spreading=np.array([1.]),
        seed=1, phase_generator="matlab",
    )
    wave = synthesize_irregular_response(
        h5_file, sea, dt=.01, end_time=10, ramp_time=20,
        rho=1025, g=9.80665,
    )
    platform = MostPlatformHydrodynamics.from_volturnus(
        h5_file, os.environ["WEC_SIM_MOST_MASS_PROPERTIES"],
    )
    tower = MostTowerReaction.from_iea15mw(
        os.environ["WEC_SIM_MOST_PROPERTIES"],
        platform_cg=platform.equilibrium_pose[:3],
    )
    # MATLAB rotor, controller, and blade-root histories remain external.
    # Tower load and platform acceleration are calculated in the Python solve.
    response = platform.simulate_with_turbine(
        time, wave.excitation_force, tower,
        source["rotor_speed"].ravel()*2*np.pi/60,
        source["azimuth"].ravel(), source["generator_torque"].ravel(),
        source["blade_aero_load"],
    )
    for axis, pose_gate, speed_gate in (
        (0, 2e-4, 5e-5),    # surge, m and m/s
        (1, 5e-4, 2e-4),    # sway, m and m/s
        (2, 1.5e-4, 7e-5),  # heave, m and m/s
        (3, 4e-5, 2e-5),    # roll, rad and rad/s
        (4, 3e-6, 2e-6),    # pitch, rad and rad/s
        (5, 1.5e-4, 2e-6),  # yaw, rad and rad/s
    ):
        np.testing.assert_allclose(
            response.position[:, axis], source["body_position"][:, axis],
            rtol=0, atol=pose_gate,
        )
        np.testing.assert_allclose(
            response.velocity[:, axis], source["body_velocity"][:, axis],
            rtol=0, atol=speed_gate,
        )

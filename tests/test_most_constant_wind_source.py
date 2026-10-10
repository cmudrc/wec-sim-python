"""Pair MOST's above-rated constant-wind platform and turbine dynamics."""

import os
from pathlib import Path

import numpy as np
import pytest
from scipy.io import loadmat

from wecsim import (
    MostBaselineController, MostConstantWind, MostCoupled,
    MostPlatformHydrodynamics, MostRotor, MostTowerReaction,
)
from wecsim.irregularWave import (
    jonswap_equal_energy_components, synthesize_irregular_response,
)

WIND_SPEED = float(os.environ.get("WEC_SIM_MOST_CONSTANT_WIND_SPEED", 12))
END_TIME = float(os.environ.get("WEC_SIM_MOST_CONSTANT_END_TIME", 10))


@pytest.mark.skipif(
    not all(os.environ.get(name) for name in (
        "WEC_SIM_MOST_CONSTANT_BASELINE", "WEC_SIM_MOST_H5",
        "WEC_SIM_MOST_MASS_PROPERTIES", "WEC_SIM_MOST_PROPERTIES",
        "WEC_SIM_MOST_BLADE_DIR", "WEC_SIM_MOST_CONTROL",
        "WEC_SIM_MOST_STEADY_STATES",
    )),
    reason="pinned MATLAB MOST constant-wind inputs not provided",
)
def test_most_above_rated_constant_wind_against_pinned_source():
    h5_file = os.environ["WEC_SIM_MOST_H5"]
    platform = MostPlatformHydrodynamics.from_volturnus(
        h5_file, os.environ["WEC_SIM_MOST_MASS_PROPERTIES"],
    )
    tower = MostTowerReaction.from_iea15mw(
        os.environ["WEC_SIM_MOST_PROPERTIES"],
        platform_cg=platform.equilibrium_pose[:3],
    )
    controller = MostBaselineController.from_matlab_files(
        os.environ["WEC_SIM_MOST_CONTROL"],
        os.environ["WEC_SIM_MOST_STEADY_STATES"], wind_speed=WIND_SPEED,
    )
    rotor = MostRotor.from_iea15mw(
        os.environ["WEC_SIM_MOST_BLADE_DIR"], controller=controller,
    )
    sea = jonswap_equal_energy_components(
        h5_file, significant_height=4, peak_period=8,
        directions=np.array([0.]), spreading=np.array([1.]),
        seed=1, phase_generator="matlab",
    )
    wave = synthesize_irregular_response(
        h5_file, sea, dt=.01, end_time=END_TIME, ramp_time=20,
        rho=1025, g=9.80665,
    )
    time = wave.time
    result = None
    if END_TIME <= 10:
        result = MostCoupled(platform, rotor, tower).simulate(
            time, wave.excitation_force,
            MostConstantWind(WIND_SPEED, end_time=END_TIME),
        )

    source = loadmat(Path(os.environ["WEC_SIM_MOST_CONSTANT_BASELINE"]))
    assert time.shape == (round(END_TIME / .01) + 1,)
    np.testing.assert_allclose(time, source["body_time"].ravel(),
                               rtol=0, atol=1e-12)
    np.testing.assert_allclose(time, source["turbine_time"].ravel(),
                               rtol=0, atol=1e-12)
    np.testing.assert_allclose(source["wind_speed"],
                               np.tile([WIND_SPEED, 0., 0.], (len(time), 1)),
                               rtol=0, atol=1e-12)
    for actual, expected in (
        (sea.omega, source["wave_omega"].ravel()),
        (sea.spectral_amplitude, source["wave_amplitude"].ravel()),
        (sea.d_omega, source["wave_d_omega"].ravel()),
        (sea.phase, source["wave_phase"]),
    ):
        np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-11)
    np.testing.assert_allclose(wave.elevation,
                               source["wave_elevation"].ravel(),
                               rtol=0, atol=1e-12)
    np.testing.assert_allclose(wave.excitation_force,
                               source["body_force_excitation"],
                               rtol=0, atol=1e-4)

    # MOST_Lib logs BladePitch after its Radians to Degrees block. BEM and
    # the Python controller use radians before that reporting conversion.
    source_pitch = np.deg2rad(source["blade_pitch"].ravel())
    assert controller.initial_pitch > .1
    assert source_pitch[0] == 0
    if result is not None:
        assert result.rotor.blade_pitch[0] == 0
    assert np.max(source_pitch) > .015
    if END_TIME > 20:
        # Exercise the developed sea while the controller keeps pitch active.
        assert np.min(source_pitch[round(20 / .01):]) > .1
    np.testing.assert_allclose(
        controller.initial_omega*60/(2*np.pi), source["rotor_speed"][0, 0],
        rtol=0, atol=1e-12,
    )
    replay_torque, replay_pitch = controller.simulate(
        time, source["rotor_speed"].ravel()*2*np.pi/60,
    )
    np.testing.assert_allclose(replay_torque,
                               source["generator_torque"].ravel(),
                               rtol=0, atol=10)
    # The sustained case's logged pitch differs slightly from a replay on
    # output-sampled rotor speed during active feedback after 24 s.
    np.testing.assert_allclose(replay_pitch, source_pitch, rtol=0,
                               atol=3e-4 if END_TIME > 20 else 1e-6)

    if result is not None:
        assert result.iterations <= 12
        assert result.position_residual <= 1e-6
        assert result.velocity_residual <= 1e-6
    causal = MostCoupled(platform, rotor, tower).simulate_causal(
        time, wave.excitation_force,
        MostConstantWind(WIND_SPEED, end_time=END_TIME),
    )
    assert causal.iterations == 1
    assert causal.position_residual <= 1e-6
    assert causal.velocity_residual <= 1e-4

    # The longer case has explicit full-duration gates, separate from the
    # 10 s ramp gates. Both use the same-run generated controller inputs.
    position_and_speed_gates = (
        (
            (0, 2e-3, 4.5e-4),
            (1, 3.5e-3, 1.2e-3),
            (2, 4.5e-4, 1.1e-4),
            (3, 4e-4, 1.2e-4),
            (4, 1.2e-4, 1.2e-5),
            (5, 5e-4, 6e-5),
        ) if END_TIME > 20 else (
            (0, 2e-4, 5e-5),
            (1, 5e-4, 2e-4),
            (2, 1.5e-4, 7e-5),
            (3, 5e-5, 2e-5),
            (4, 5e-6, 2e-6),
            (5, 1.5e-4, 3e-6),
        )
    )
    for coupled in ((result, causal) if result is not None else (causal,)):
        for axis, position_gate, velocity_gate in position_and_speed_gates:
            np.testing.assert_allclose(
                coupled.platform.position[:, axis],
                source["body_position"][:, axis],
                rtol=0, atol=position_gate,
            )
            np.testing.assert_allclose(
                coupled.platform.velocity[:, axis],
                source["body_velocity"][:, axis],
                rtol=0, atol=velocity_gate,
            )
        np.testing.assert_allclose(
            coupled.rotor.rotor_speed*60/(2*np.pi),
            source["rotor_speed"].ravel(), rtol=0,
            atol=.004 if END_TIME > 20 else .001,
        )
        np.testing.assert_allclose(
            coupled.rotor.azimuth, source["azimuth"].ravel(),
            rtol=0, atol=.006 if END_TIME > 20 else .0003,
        )
        np.testing.assert_allclose(
            coupled.rotor.generator_torque,
            source["generator_torque"].ravel(), rtol=0,
            atol=100 if END_TIME > 20 else 3000,
        )
        np.testing.assert_allclose(
            coupled.rotor.blade_pitch, source_pitch, rtol=0,
            atol=.001 if END_TIME > 20 else 2e-5,
        )
        expected_load = source["blade_aero_load"]
        assert coupled.rotor.blade_root_load.shape == expected_load.shape == (
            len(time), 6, 3,
        )
        component_peak = np.max(np.abs(expected_load), axis=(0, 2))
        component_error = np.max(
            np.abs(coupled.rotor.blade_root_load - expected_load), axis=(0, 2),
        )
        assert np.all(component_error < (
            .006 if END_TIME > 20 else .002
        ) * component_peak), (
            component_error, component_peak,
        )

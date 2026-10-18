"""Independently advance the published WaveStar plant before NMPC starts."""

from pathlib import Path

import numpy as np
import pytest

from wecsim import WaveStarNmpcPTO
from wecsim.dynamics.wavestar import WaveStarLinkage, run_wavestar_published
from wecsim.environment.irregularWave import jonswap_equal_energy_components


ROOT = Path(__file__).parent / "fixtures/matlab"
HYDRO = ROOT / "wavestar_fixed_input_step/wavestar.h5"
SOURCE = ROOT / "wavestar_fixed_input_step/source.npz"
FULL_SOURCE = ROOT / "wavestar_nmpc_full_step/baseline.npz"
AUXILIARY = ROOT / "wavestar_auxiliary_mass/source.npz"
ORIGINAL_AUXILIARY = ROOT / "wavestar_added_mass_delay/source.npz"


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


def test_light_auxiliary_bodies_have_bounded_precontrol_effect():
    """A source run reduces BC/motor mass and inertia with the same 50 ms sea."""
    with np.load(SOURCE, allow_pickle=False) as published, \
            np.load(AUXILIARY, allow_pickle=False) as auxiliary, \
            np.load(ORIGINAL_AUXILIARY, allow_pickle=False) as original:
        sea = _published_sea(published["50ms_components"])
        source_time = published[
            "50ms_WECCCOMP_Nonlinear_Model_Predictive_body1"
        ][:, 0]
        assert auxiliary["auxiliary_mass"] == 1e-9
        np.testing.assert_allclose(auxiliary["time"], source_time,
                                   rtol=0, atol=1e-12)
        initial_change = (auxiliary["acceleration"][0, 4]
                          - original["acceleration"][0, 4])
        assert 3.8e-5 < initial_change < 4e-5
        assert np.max(np.abs(
            auxiliary["pitch"] - original["pitch"]
        )) < 3e-6
        assert np.max(np.abs(
            auxiliary["speed"] - original["speed"]
        )) < 1.1e-5
        response = run_wavestar_published(
            HYDRO, sea, dt=.05, end_time=9.95, ramp_time=25,
            g=9.80665, solver="wecsim_ode8",
        )
        # This is a diagnostic source variant. The published case is gated
        # separately above, and its original small bodies remain in source.
        assert np.max(np.abs(response.angle - auxiliary["pitch"])) < 2e-5
        assert np.max(np.abs(
            response.angular_speed - auxiliary["speed"]
        )) < 1.5e-4


def test_auxiliary_gravity_dominates_the_initial_acceleration_shift():
    """Project the published light bodies onto the initial arm coordinate."""
    with np.load(ORIGINAL_AUXILIARY, allow_pickle=False) as published, \
            np.load(AUXILIARY, allow_pickle=False) as near_zero:
        assert published["pitch"][0] == 0
        source_acceleration = published["acceleration"][0, 4]
        near_zero_acceleration = near_zero["acceleration"][0, 4]

    linkage = WaveStarLinkage()
    base = np.asarray(linkage.pto_base)
    tip = np.asarray(linkage.pto_arm)
    rod = tip - base
    tip_tangent = linkage.point_velocity(tip, 0., 1.)
    rod_angle_rate = (rod[1] * tip_tangent[0] - rod[0] * tip_tangent[1]) \
        / (rod @ rod)
    rod_tangent = tip_tangent - rod_angle_rate * np.array([tip[1], -tip[0]])
    motor_tangent = -rod_angle_rate * np.array([base[1], -base[0]])
    projected_mass = 1e-4 * (
        rod_tangent @ rod_tangent + motor_tangent @ motor_tangent
        + 2 * rod_angle_rate**2
    )
    gravity_load = 1e-4 * 9.80665 * (rod_tangent[1] + motor_tangent[1])

    # The gravity load is over 100 times the initial inertia force. Both
    # together imply the same approximately 3.33 kg m² base mass as the
    # independent float/arm hydrodynamic model.
    assert gravity_load > 100 * projected_mass * near_zero_acceleration
    base_mass = (gravity_load + source_acceleration * projected_mass) \
        / (near_zero_acceleration - source_acceleration)
    assert 3.32 < base_mass < 3.34


def test_source_ode8_precontrol_motion_on_two_refined_steps():
    with np.load(SOURCE, allow_pickle=False) as published:
        sea = _published_sea(published["50ms_components"])
        source_bodies = {
            label: published[
                f"{label}_WECCCOMP_Nonlinear_Model_Predictive_body1"
            ] for label in ("25ms", "12p5ms")
        }
    errors = {}
    for label, dt, stride, pitch_limit, speed_limit in (
        ("25ms", .025, 2, 1e-5, 3e-5),
        ("12p5ms", .0125, 4, 5e-6, 1.5e-5),
    ):
        response = run_wavestar_published(
            HYDRO, sea, dt=dt, end_time=9.95, ramp_time=25,
            g=9.80665, solver="wecsim_ode8", output_stride=stride,
        )
        body = source_bodies[label]
        np.testing.assert_allclose(response.time, body[:, 0],
                                   rtol=0, atol=1e-10)
        pitch_gap = np.max(np.abs(response.angle - body[:, 5]))
        speed_gap = np.max(np.abs(response.angular_speed - body[:, 11]))
        assert pitch_gap < pitch_limit
        assert speed_gap < speed_limit
        errors[label] = pitch_gap, speed_gap
    assert np.all(np.array(errors["12p5ms"]) < np.array(errors["25ms"]))


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

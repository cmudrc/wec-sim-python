"""Check moving Morison force and mass behavior without source artifacts."""

import numpy as np

from wecsim.hydrodynamics.morison import (
    MorisonElement,
    regular_wave_axial_morison_terms,
    regular_wave_normal_tangential_morison_terms,
)


def _state(*, z=-2., surge=0., heave=0., pitch=0., yaw=0.):
    return np.array([surge, 0., z + heave, 0., pitch, yaw])


def test_axial_still_water_drag_and_added_mass_follow_body_velocity():
    element = MorisonElement(
        point=(0, 0, -1), drag_coefficient=(0, 0, 2),
        added_mass_coefficient=(0, 0, 1.5), area=(0, 0, 3), volume=2,
    )
    wrench, mass = regular_wave_axial_morison_terms(
        [element], position=_state(), velocity=np.array([0, 0, 1, 0, 0, 0]),
        time=0, wave_height=0, wave_period=8, ramp_time=0,
        water_depth=np.inf, rho=1000,
    )
    np.testing.assert_allclose(wrench, [0, 0, -3000, 0, 0, 0], atol=1e-12)
    np.testing.assert_allclose(mass[2, 2], 3000)
    np.testing.assert_allclose(mass, mass.T)
    assert np.linalg.eigvalsh(mass).min() >= -1e-10

    # A dry element contributes neither drag nor added mass.
    dry_wrench, dry_mass = regular_wave_axial_morison_terms(
        [element], position=_state(z=2), velocity=np.array([0, 0, 1, 0, 0, 0]),
        time=0, wave_height=0, wave_period=8, ramp_time=0,
        water_depth=np.inf, rho=1000,
    )
    np.testing.assert_array_equal(dry_wrench, 0)
    np.testing.assert_array_equal(dry_mass, 0)


def test_normal_tangential_drag_rotates_with_body_and_mass_stays_passive():
    element = MorisonElement(
        point=(0, 0, -1), axis=(1, 0, 0),
        drag_coefficient=(2, 3, 0),
        added_mass_coefficient=(1, 1.5, 0),
        area=(4, 5, 0), volume=2,
    )
    velocity = np.array([1., 0., 1., 0., 0., 0.])
    wrench, mass = regular_wave_normal_tangential_morison_terms(
        [element], position=_state(), velocity=velocity,
        time=0, wave_height=0, wave_period=8, ramp_time=0,
        water_depth=np.inf, rho=1000,
    )
    np.testing.assert_allclose(wrench[:3], [-7500, 0, -4000], atol=1e-12)
    np.testing.assert_allclose(mass, mass.T, atol=1e-12)
    assert np.linalg.eigvalsh(mass).min() >= -1e-8

    rotated, _ = regular_wave_normal_tangential_morison_terms(
        [element], position=_state(yaw=np.pi / 2),
        velocity=np.array([0., 1., 1., 0., 0., 0.]),
        time=0, wave_height=0, wave_period=8, ramp_time=0,
        water_depth=np.inf, rho=1000,
    )
    np.testing.assert_allclose(rotated[:3], [0, -7500, -4000], atol=1e-12)


def test_finite_depth_wave_and_current_yield_finite_moving_wrench():
    element = MorisonElement(
        point=(0, 0, -1), axis=(1, 0, 0),
        drag_coefficient=(2, 3, 0),
        added_mass_coefficient=(1, 1.5, 0),
        area=(4, 5, 0), volume=2,
    )
    wrench, mass = regular_wave_normal_tangential_morison_terms(
        [element], position=_state(surge=.2, pitch=.1),
        velocity=np.array([.1, 0, .2, 0, .05, 0]),
        time=1, wave_height=1, wave_period=8, direction=30,
        ramp_time=2, water_depth=10, rho=1000,
        current_speed=.5, current_direction=90,
        current_profile="linear", current_depth=5,
    )
    assert np.isfinite(wrench).all()
    assert np.linalg.norm(wrench) > 0
    assert np.linalg.eigvalsh(mass).min() >= -1e-8

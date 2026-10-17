"""Check moving Morison force and mass behavior without source artifacts."""

import numpy as np

from wecsim.hydrodynamics.morison import (
    MorisonElement,
    no_wave_heave_morison_terms,
    regular_wave_axial_morison_terms,
    regular_wave_heave_morison_terms,
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


def test_pitched_axial_element_keeps_rigid_geometry_and_resists_relative_motion():
    """An axial element rotates without changing length or creating a spurious moment."""
    pitch = .3
    axis = np.array([np.sin(pitch), 0., np.cos(pitch)])
    element = MorisonElement(
        point=(0, 0, -1), drag_coefficient=(0, 0, 2),
        added_mass_coefficient=(0, 0, 0), area=(0, 0, 3), volume=1,
    )
    settings = dict(
        position=_state(pitch=pitch), time=1, wave_height=0,
        wave_period=8, ramp_time=0, water_depth=100, rho=1000,
    )
    current, current_mass = regular_wave_axial_morison_terms(
        [element], velocity=np.zeros(6), current_speed=1, **settings,
    )
    # The body-local one-meter arm and drag force both follow the pitched
    # axis. Horizontal current has axial speed sin(pitch), so no moment is
    # generated about the body center for this collinear element.
    expected = 3000 * np.sin(pitch)**2 * axis
    np.testing.assert_allclose(current[:3], expected, rtol=0, atol=1e-10)
    np.testing.assert_allclose(current[3:], 0, rtol=0, atol=1e-10)
    np.testing.assert_array_equal(current_mass, 0)

    body_velocity = np.array([.8, 0, 0, 0, 0, 0])
    still, _ = regular_wave_axial_morison_terms(
        [element], velocity=body_velocity, **settings,
    )
    np.testing.assert_allclose(
        still[:3], -3000 * (.8 * np.sin(pitch))**2 * axis,
        rtol=0, atol=1e-10,
    )
    np.testing.assert_allclose(still[3:], 0, rtol=0, atol=1e-10)
    assert still[:3] @ body_velocity[:3] < 0


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


def test_heave_element_dries_at_mean_waterline_and_bottom_has_no_vertical_orbit():
    element = MorisonElement(
        point=(0, 0, 0), drag_coefficient=(0, 0, 2),
        added_mass_coefficient=(0, 0, 1), area=(0, 0, 3), volume=2,
    )
    # A submerged element resists upward motion; just above the mean
    # waterline both drag and added mass switch off.
    assert no_wave_heave_morison_terms(
        [element], center_z=0, heave=0, speed=1, rho=1000,
    ) == (-3000, 2000)
    assert no_wave_heave_morison_terms(
        [element], center_z=.01, heave=0, speed=1, rho=1000,
    ) == (0, 0)

    settings = dict(time=0, wave_height=2, wave_period=8, ramp_time=0,
                    water_depth=10, rho=1000)
    # Linear-wave vertical orbital motion vanishes at the seabed. The
    # submerged element still contributes its acceleration-dependent mass.
    assert regular_wave_heave_morison_terms(
        [element], center_z=-10, heave=0, speed=0, **settings,
    ) == (0, 2000)
    mid_force, mid_mass = regular_wave_heave_morison_terms(
        [element], center_z=-5, heave=0, speed=0, **settings,
    )
    assert mid_force < 0 and mid_mass == 2000
    assert regular_wave_heave_morison_terms(
        [element], center_z=.01, heave=0, speed=0, **settings,
    ) == (0, 0)

    wrench, mass = regular_wave_axial_morison_terms(
        [element], position=(0, 0, -10, 0, 0, 0),
        velocity=np.zeros(6), **settings,
    )
    np.testing.assert_allclose(wrench, 0, atol=1e-12)
    assert mass[2, 2] == 2000


def test_normal_tangential_element_dries_with_body_motion():
    element = MorisonElement(
        point=(0, 0, -1), axis=(1, 0, 0),
        drag_coefficient=(2, 3, 0),
        added_mass_coefficient=(1, 1.5, 0),
        area=(4, 5, 0), volume=2,
    )
    settings = dict(velocity=(1, 0, 1, 0, 0, 0), time=0,
                    wave_height=0, wave_period=8, ramp_time=0,
                    water_depth=np.inf, rho=1000)
    wet_wrench, wet_mass = regular_wave_normal_tangential_morison_terms(
        [element], position=(0, 0, 0, 0, 0, 0), **settings,
    )
    assert np.linalg.norm(wet_wrench) > 0
    assert np.trace(wet_mass) > 0
    dry_wrench, dry_mass = regular_wave_normal_tangential_morison_terms(
        [element], position=(0, 0, 1.01, 0, 0, 0), **settings,
    )
    np.testing.assert_array_equal(dry_wrench, 0)
    np.testing.assert_array_equal(dry_mass, 0)

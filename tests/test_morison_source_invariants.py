"""Known simple limits of the MATLAB Morison source-law diagnostics."""

import numpy as np

from wecsim.environment.irregularWave import IrregularComponents
from wecsim.hydrodynamics.morison import (
    MorisonElement, finite_depth_wavenumber,
    irregular_morison_source_drag, regular_morison_source_force,
)


def _state(z=-2.):
    return np.array([0., 0., z, 0., 0., 0.])


def test_regular_cartesian_source_law_reduces_to_drag_plus_added_mass():
    element = MorisonElement(
        point=(0, 0, -1), drag_coefficient=(0, 0, 2),
        added_mass_coefficient=(0, 0, 1.5), area=(0, 0, 3), volume=2,
    )
    force = regular_morison_source_force(
        [element], time=0, position=_state(),
        velocity=np.array([0, 0, 1, 0, 0, 0]),
        acceleration=np.array([0, 0, .2, 0, 0, 0]),
        wave_height=0, wave_period=8, direction=0,
        water_depth=10, ramp_time=0, rho=1000,
    )
    np.testing.assert_allclose(force, [0, 0, -3600, 0, 0, 0], atol=1e-10)

    dry = regular_morison_source_force(
        [element], time=0, position=_state(2),
        velocity=np.array([0, 0, 1, 0, 0, 0]),
        acceleration=np.zeros(6), wave_height=0, wave_period=8,
        direction=0, water_depth=10, ramp_time=0, rho=1000,
    )
    np.testing.assert_array_equal(dry, 0)


def test_regular_normal_tangential_source_law_splits_drag_components():
    element = MorisonElement(
        point=(0, 0, -1), drag_coefficient=(2, 3, 0),
        added_mass_coefficient=(0, 0, 0), area=(4, 5, 0), volume=2,
    )
    force = regular_morison_source_force(
        [element], time=0, position=_state(),
        velocity=np.array([1, 0, 1, 0, 0, 0]), acceleration=np.zeros(6),
        wave_height=0, wave_period=8, direction=0,
        water_depth=10, ramp_time=0, rho=1000,
        body_morison=2, element_axes=[(1, 0, 0)],
    )
    np.testing.assert_allclose(force[:3], [-7500, 0, -4000], atol=1e-10)
    np.testing.assert_allclose(force[4], 7500, atol=1e-10)


def test_irregular_source_drag_matches_linear_wave_velocity_at_rest():
    omega = np.array([.8, 1.2])
    spectrum = np.array([.2, .1])
    width = np.array([.1, .1])
    components = IrregularComponents(
        omega=omega, spectral_amplitude=spectrum, d_omega=width,
        directions=np.array([0.]), spreading=np.array([1.]),
        phase=np.zeros((2, 1)),
    )
    element = MorisonElement(
        point=(0, 0, -1), drag_coefficient=(2, 0, 0),
        added_mass_coefficient=(0, 0, 0), area=(3, 0, 0), volume=0,
    )
    force = irregular_morison_source_drag(
        [element], time=0, position=_state(), velocity=np.zeros(6),
        components=components, water_depth=10, ramp_time=0, rho=1000,
    )
    k = finite_depth_wavenumber(omega, water_depth=10, gravity=9.81)
    horizontal = np.cosh(k * 7) / np.cosh(k * 10)
    speed = np.sum(np.sqrt(spectrum * width) * horizontal * 9.81 * k / omega)
    np.testing.assert_allclose(force[:3], [3000 * speed**2, 0, 0],
                               rtol=0, atol=1e-10)
    np.testing.assert_allclose(force[4], -3000 * speed**2, rtol=0, atol=1e-10)

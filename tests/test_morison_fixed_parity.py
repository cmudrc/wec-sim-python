"""Pair the published no-HDF5 Morison monopile with pinned MATLAB output."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim import Current, PMWave, WEC
from wecsim.environment.irregularWave import IrregularComponents, pm_equal_energy_components
from wecsim.hydrodynamics.morison import (
    MorisonElement, finite_depth_wavenumber, solve_fixed_morison_irregular,
)


REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")


def _max_error(actual, expected, limit, label):
    actual = np.asarray(actual)
    expected = np.asarray(expected)
    assert actual.shape == expected.shape, label
    assert np.isfinite(actual).all() and np.isfinite(expected).all(), label
    error = np.max(np.abs(actual - expected))
    assert error < limit, f"{label}: {error:.6g} exceeds {limit}"


@pytest.mark.skipif(
    not REFERENCE, reason="paired MATLAB fixed Morison output not provided",
)
def test_published_fixed_morison_wave_and_force():
    source = Path(REFERENCE)
    components_csv = np.loadtxt(source / "MORISON_FIXED_components.csv", delimiter=",")
    directions = np.loadtxt(source / "MORISON_FIXED_directions.csv", delimiter=",")
    wave = np.loadtxt(source / "MORISON_FIXED_wave.csv", delimiter=",")
    force = np.loadtxt(source / "MORISON_FIXED_morison_force.csv", delimiter=",")
    monopile = np.loadtxt(source / "MORISON_FIXED_morisonElement_body1.csv",
                          delimiter=",")
    tower = np.loadtxt(source / "MORISON_FIXED_morisonElement_body2.csv",
                       delimiter=",")
    assert components_csv.shape == (500, 7)
    assert directions.shape == (3, 2)
    assert wave.shape == (40001, 2)
    assert force.shape == (40001, 7)
    assert monopile.shape == tower.shape == (40001, 25)
    _max_error(monopile[:, 13:19], -force[:, 1:7], 1e-8,
               "source total force and signed Morison output")

    components = pm_equal_energy_components(
        None, significant_height=2, peak_period=5,
        directions=directions[:, 0], spreading=directions[:, 1],
        frequency_range=(.001, 10), seed=5, phase_generator="matlab",
    )
    _max_error(components.phase, components_csv[:, 4:7], 1e-13,
               "seeded three-direction MATLAB phases")
    for actual, expected, label in (
        (components.omega, components_csv[:, 0], "frequency"),
        (components.spectral_amplitude, components_csv[:, 1], "PM spectrum"),
        (components.d_omega, components_csv[:, 2], "frequency width"),
        (finite_depth_wavenumber(components.omega, water_depth=30),
         components_csv[:, 3], "finite-depth wavenumber"),
    ):
        _max_error(actual, expected, 1e-12, label)

    wec = WEC("published fixed Morison monopile")
    body = wec.fixed_body(
        "monopile", center_gravity=(0, 0, -15),
        mass="equilibrium", inertia=(1.25e9, 1.25e9, .15e9),
        volume=np.pi * 10**2 * 30,
    )
    wec.fixed_body(
        "tower", center_gravity=(0, 0, 25), mass=1031930,
        inertia=(9.66e8, 9.66e8, .132e8),
    )
    wec.morison_element(
        body, point=body.at(0, 0, 10),
        drag_coefficient=(1, 1, 1),
        added_mass_coefficient=(1, 1, 1),
        area=(300, 300, np.pi * 10**2 / 4),
        volume=np.pi * 10**2 * 30,
        phase_mode="matlab_shared",
    )
    result = wec.run(
        PMWave(2, 5, seed=5, phase_generator="matlab",
               directions=(0, 30, 90), spreading=(.1, .2, .7),
               frequency_range=(.001, 10), water_depth=30),
        dt=.01, end_time=400, ramp_time=100, rho=1025,
    )
    _max_error(result.time, wave[:, 0], 1e-12, "time")
    _max_error(result.wave_elevation, wave[:, 1], 1e-10,
               "published wave elevation")
    _max_error(result.body_forces["monopile"], -force[:, 1:7], 1e-3,
               "published six-component Morison force")
    _max_error(result.body_forces["tower"], np.zeros((40001, 6)), 1e-12,
               "tower Morison force")
    for name, saved in (("monopile", monopile), ("tower", tower)):
        _max_error(result.bodies[name].position, saved[:, 1:7], 1e-12,
                   f"{name} position")
        _max_error(result.bodies[name].velocity, saved[:, 7:13], 1e-12,
                   f"{name} velocity")


def test_directional_morison_phase_is_the_default():
    components = IrregularComponents(
        omega=np.array([1., 2.]),
        spectral_amplitude=np.array([.5, .3]),
        d_omega=np.array([.1, .2]),
        directions=np.array([0., 90.]),
        spreading=np.array([.5, .5]),
        phase=np.array([[0., np.pi / 2], [0., np.pi / 2]]),
    )
    common = dict(
        point=(0, 0, 0), drag_coefficient=(0, 0, 0),
        added_mass_coefficient=(0, 0, 0), area=(0, 0, 0), volume=1,
    )
    options = dict(center_gravity=(0, 0, -5), water_depth=30,
                   dt=.1, end_time=.1, ramp_time=0, rho=1)
    directional = solve_fixed_morison_irregular(
        components, [MorisonElement(**common)], **options,
    )
    source_shared = solve_fixed_morison_irregular(
        components, [MorisonElement(**common, phase_mode="matlab_shared")],
        **options,
    )
    _max_error(directional.wave_elevation, source_shared.wave_elevation,
               1e-12, "phase mode must not alter the sea")
    assert np.max(np.abs(directional.force - source_shared.force)) > .01


def test_fixed_body_default_zero_inertia_matches_explicit_stationary_inertia():
    """A stationary Morison body has no inertia term in its force equation."""
    def run(inertia=None):
        wec = WEC("fixed pile")
        settings = {} if inertia is None else {"inertia": inertia}
        pile = wec.fixed_body("pile", center_gravity=(0, 0, -5),
                              volume=1, **settings)
        wec.morison_element(
            pile, point=pile.at(0, 0, 0),
            drag_coefficient=(1, 1, 1),
            added_mass_coefficient=(0, 0, 0),
            area=(2, 2, 2), volume=1,
        )
        return wec.run(
            PMWave(2, 5, seed=1, frequency_count=4,
                   frequency_range=(.1, 2), water_depth=30),
            dt=.1, end_time=.2, ramp_time=0,
        )

    default = run()
    explicit = run((1, 2, 3))
    assert default.case["bodies"][0]["inertia"] == [0, 0, 0]
    assert np.max(np.abs(default.body_forces["pile"])) > 0
    np.testing.assert_array_equal(default.body_forces["pile"],
                                  explicit.body_forces["pile"])
    np.testing.assert_array_equal(default.bodies["pile"].position,
                                  explicit.bodies["pile"].position)
    np.testing.assert_array_equal(default.wave_elevation,
                                  explicit.wave_elevation)
    with pytest.raises(ValueError, match="three nonnegative values"):
        run((-1, 0, 0))


@pytest.mark.parametrize("profile,depth,scale", [
    ("uniform", None, 1),
    ("power", 10, .5**(1 / 7)),
    ("linear", 10, .5),
    ("linear", 4, 0),
])
def test_fixed_morison_current_drag_profile_and_ramp(profile, depth, scale):
    components = IrregularComponents(
        omega=np.array([1., 2.]), spectral_amplitude=np.zeros(2),
        d_omega=np.ones(2), directions=np.array([0.]),
        spreading=np.ones(1), phase=np.zeros((2, 1)),
    )
    element = MorisonElement(
        point=(0, 0, 10), drag_coefficient=(1, 0, 0),
        added_mass_coefficient=(0, 0, 0), area=(2, 0, 0), volume=1,
    )
    result = solve_fixed_morison_irregular(
        components, [element], center_gravity=(0, 0, -15),
        water_depth=30, dt=1, end_time=2, ramp_time=2, rho=1000,
        current_speed=2, current_profile=profile, current_depth=depth,
    )
    full_force = .5 * 1000 * 2 * (2 * scale)**2
    np.testing.assert_allclose(result.force[:, 0],
                               full_force * np.array([0, .25, 1]), atol=1e-12)
    np.testing.assert_allclose(result.force[:, 4], 10 * result.force[:, 0])
    np.testing.assert_allclose(result.force[:, [1, 2, 3, 5]], 0)


def test_current_public_api_and_multiheading_limit():
    wave = PMWave(2, 5, current=Current(.8, 45, "power", 30))
    assert wave.as_case()["current"] == {
        "speed": .8, "direction": 45, "profile": "power", "depth": 30,
    }
    with pytest.raises(ValueError, match="depth"):
        Current(.8, profile="power").as_case()
    components = IrregularComponents(
        omega=np.array([1., 2.]), spectral_amplitude=np.zeros(2),
        d_omega=np.ones(2), directions=np.array([0., 90.]),
        spreading=np.array([.5, .5]), phase=np.zeros((2, 2)),
    )
    element = MorisonElement(
        point=(0, 0, 0), drag_coefficient=(1, 0, 0),
        added_mass_coefficient=(0, 0, 0), area=(1, 0, 0), volume=1,
    )
    with pytest.raises(ValueError, match="one incident heading"):
        solve_fixed_morison_irregular(
            components, [element], center_gravity=(0, 0, -5),
            water_depth=30, dt=.1, end_time=.1, ramp_time=0,
            current_speed=.8,
        )

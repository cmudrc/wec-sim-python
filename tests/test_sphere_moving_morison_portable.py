"""Portable checks for the moving-Morison Sphere path paired in MATLAB CI."""

from pathlib import Path

import numpy as np

from wecsim import NoWave, WEC
from wecsim.hydrodynamics.morison import (
    MorisonElement,
    no_wave_heave_morison_terms,
    regular_wave_axial_morison_terms,
)


SPHERE = Path(__file__).resolve().parents[1] / "examples/data/sphere.h5"
MATLAB = (Path(__file__).resolve().parent / "fixtures/matlab"
          / "sphere_moving_morison_heave.npz")


def _element(drag_coefficient=1):
    return MorisonElement(
        point=(0, 0, -2), drag_coefficient=(0, 0, drag_coefficient),
        added_mass_coefficient=(0, 0, 1), area=(0, 0, 100), volume=20,
    )


def test_submerged_axial_element_matches_still_water_force_and_mass():
    element = _element()
    drag, added = no_wave_heave_morison_terms(
        [element], center_z=-2, heave=1, speed=1, rho=1000,
    )
    assert drag == -50_000
    assert added == 20_000

    wrench, mass = regular_wave_axial_morison_terms(
        [element], position=(0, 0, -1, 0, 0, 0),
        velocity=(0, 0, 1, 0, 0, 0), time=0,
        wave_height=0, wave_period=8, ramp_time=0,
        water_depth=np.inf, rho=1000,
    )
    np.testing.assert_allclose(wrench, (0, 0, drag, 0, 0, 0), atol=1e-12)
    expected_mass = np.zeros((6, 6))
    expected_mass[2, 2] = added
    np.testing.assert_allclose(mass, expected_mass, atol=1e-12)

    # Lifting the attachment point clear of the water removes both terms.
    dry_wrench, dry_mass = regular_wave_axial_morison_terms(
        [element], position=(0, 0, 3, 0, 0, 0),
        velocity=(0, 0, 1, 0, 0, 0), time=0,
        wave_height=0, wave_period=8, ramp_time=0,
        water_depth=np.inf, rho=1000,
    )
    np.testing.assert_array_equal(dry_wrench, 0)
    np.testing.assert_array_equal(dry_mass, 0)


def _run_sphere(drag_coefficient):
    wec = WEC("moving Morison sphere")
    sphere = wec.body("sphere", SPHERE)
    wec.coordinate("heave", sphere.move("heave"))
    wec.morison_element(
        sphere, point=sphere.at(0, 0, -2),
        drag_coefficient=(0, 0, drag_coefficient),
        added_mass_coefficient=(0, 0, 1), area=(0, 0, 100), volume=20,
    )
    return wec.run(
        NoWave(), dt=.01, end_time=.1, radiation_memory=.1,
        initial_coordinate={"heave": 1}, initial_speed={"heave": 1},
    )


def test_public_sphere_run_applies_drag_to_its_own_motion():
    undamped = _run_sphere(0)
    damped = _run_sphere(1)
    np.testing.assert_allclose(damped.time, np.arange(11) * .01, atol=1e-12)
    assert damped.wave_elevation is None
    np.testing.assert_allclose(
        damped.bodies["sphere"].position[0, 2], -1, atol=1e-12,
    )
    np.testing.assert_allclose(
        damped.bodies["sphere"].velocity[0, 2], 1, atol=1e-12,
    )
    # Both runs have identical added mass. Drag must oppose the initial
    # upward speed and reduce it over the short, otherwise identical run.
    assert (damped.body_forces["sphere"][0, 2]
            < undamped.body_forces["sphere"][0, 2] - 40_000)
    assert (damped.bodies["sphere"].velocity[-1, 2]
            < undamped.bodies["sphere"].velocity[-1, 2] - .005)
    assert (damped.bodies["sphere"].position[-1, 2]
            < undamped.bodies["sphere"].position[-1, 2])
    assert np.isfinite(damped.body_forces["sphere"]).all()


def test_full_duration_sphere_heave_against_pinned_matlab():
    """Independently advance the 40 s derived source case in ordinary CI."""
    wec = WEC("Sphere with moving axial Morison element")
    sphere = wec.body("sphere", SPHERE)
    wec.coordinate("heave", sphere.move("heave"))
    wec.morison_element(
        sphere, point=sphere.at(0, 0, -2),
        drag_coefficient=(0, 0, 1), added_mass_coefficient=(0, 0, 1),
        area=(0, 0, 100), volume=20,
    )
    result = wec.run(
        NoWave(), dt=.01, end_time=40, radiation_memory=15,
        initial_coordinate={"heave": 1},
    )
    with np.load(MATLAB, allow_pickle=False) as source:
        assert source["time"].shape == (4001,)
        np.testing.assert_allclose(result.time, source["time"],
                                   rtol=0, atol=1e-12)
        heave_error = np.max(np.abs(
            result.bodies["sphere"].position[:, 2]
            - source["position_heave"],
        ))
        speed_error = np.max(np.abs(
            result.bodies["sphere"].velocity[:, 2]
            - source["velocity_heave"],
        ))
        force_error = np.max(np.abs(
            result.body_forces["sphere"][100:, 2]
            - source["physical_force_heave"][100:],
        ))
        force_impulse_error = abs(np.trapezoid(
            result.body_forces["sphere"][:, 2]
            - source["physical_force_heave"], result.time,
        ))
        assert heave_error < .0005, heave_error
        assert speed_error < .0011, speed_error
        assert force_error < 50, force_error
        assert force_impulse_error < 250, force_impulse_error
        # The source logs zero acceleration-feedback force at startup;
        # this known transient is excluded only from the pointwise gate.
        assert abs(result.body_forces["sphere"][0, 2]
                   - source["physical_force_heave"][0]) > 30_000

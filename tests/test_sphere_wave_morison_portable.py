"""Full-duration moving Sphere Morison gate from pinned R2025b output."""

from pathlib import Path

import numpy as np

from wecsim import RegularWave, WEC
from wecsim.hydrodynamics.morison import (
    MorisonElement, regular_wave_heave_morison_terms,
)


FIXTURE = Path(__file__).parent / "fixtures/matlab/sphere_wave_morison"


def _source():
    with np.load(FIXTURE / "source.npz", allow_pickle=False) as saved:
        body, wave, force = (saved[name] for name in ("body", "wave", "source_force"))
    assert body.shape == (4001, 37)
    assert wave.shape == (4001, 2)
    assert force.shape == (4001, 3)
    assert all(np.isfinite(values).all() for values in (body, wave, force))
    time = np.arange(4001) * .01
    for recorded in (body[:, 0], wave[:, 0], force[:, 0]):
        np.testing.assert_allclose(recorded, time, rtol=0, atol=1e-12)
    return body, wave, force


def _within(actual, expected, limit, label):
    actual = np.asarray(actual)
    expected = np.asarray(expected)
    assert actual.shape == expected.shape, label
    difference = float(np.max(np.abs(actual - expected)))
    assert np.isfinite(difference) and difference < limit, (
        f"{label}: {difference:.6g} exceeds {limit:.6g}"
    )


def test_published_morison_law_on_projected_source_heave():
    body, _, source_force = _source()
    element = MorisonElement(
        point=(0, 0, -2), drag_coefficient=(0, 0, 1),
        added_mass_coefficient=(0, 0, 1), area=(0, 0, 100), volume=20,
    )
    applied_and_mass = np.array([
        regular_wave_heave_morison_terms(
            [element], center_z=-2, heave=row[3] + 2,
            speed=row[9], time=row[0], wave_height=1, wave_period=8,
            ramp_time=10, water_depth=np.inf, rho=1000,
        )
        for row in body
    ])
    reconstructed = applied_and_mass[:, 0] - applied_and_mass[:, 1] * body[:, 33]
    _within(reconstructed, source_force[:, 2], 1e-6,
            "Morison law on the source's projected heave")
    _within(source_force[50:, 1], -body[50:, 27], 20,
            "direct source force after startup")
    assert np.max(np.abs(body[:, 1])) > 1


def test_independent_public_sphere_heave_and_force_over_40_seconds():
    body, wave, _ = _source()
    wec = WEC("Sphere moving Morison in regular waves")
    sphere = wec.body("sphere", FIXTURE / "sphere.h5")
    wec.coordinate("heave", sphere.move("heave"))
    wec.morison_element(
        sphere, point=sphere.at(0, 0, -2),
        drag_coefficient=(0, 0, 1),
        added_mass_coefficient=(0, 0, 1),
        area=(0, 0, 100), volume=20,
    )
    result = wec.run(
        RegularWave(1, 8), dt=.01, end_time=40, ramp_time=10,
        initial_coordinate={"heave": 1},
    )
    _within(result.time, body[:, 0], 1e-12, "time")
    _within(result.wave_elevation, wave[:, 1], 1e-10, "wave elevation")
    _within(result.bodies["sphere"].position[:, 2], body[:, 3], .01,
            "heave position")
    _within(result.bodies["sphere"].velocity[:, 2], body[:, 9], .02,
            "heave speed")
    _within(result.body_forces["sphere"][100:, 2], -body[100:, 27], 1000,
            "moving Morison force after startup")


def test_independent_public_sphere_three_coordinate_motion_and_force():
    """Gate the published planar joint without prescribing source motion."""
    body, wave, _ = _source()
    wec = WEC("Sphere three-coordinate moving Morison")
    sphere = wec.body(
        "sphere", FIXTURE / "sphere.h5",
        inertia=(20907301, 21306090.66, 37085481.11),
    )
    for axis in ("surge", "heave", "pitch"):
        wec.coordinate(axis, sphere.move(axis))
    wec.morison_element(
        sphere, point=sphere.at(0, 0, -2),
        drag_coefficient=(0, 0, 1),
        added_mass_coefficient=(0, 0, 1),
        area=(0, 0, 100), volume=20,
    )
    result = wec.run(
        RegularWave(1, 8), dt=.01, end_time=40, ramp_time=10,
        initial_coordinate={"heave": 1},
    )
    _within(result.time, body[:, 0], 1e-12, "time")
    _within(result.wave_elevation, wave[:, 1], 1e-10, "wave elevation")
    for axis, name, position_limit, speed_limit, force_limit in (
        (0, "surge", .015, .002, 500),
        (2, "heave", .0015, .002, 200),
        (4, "pitch", .00015, .00006, 250),
    ):
        _within(result.bodies["sphere"].position[:, axis],
                body[:, 1 + axis], position_limit, f"{name} position")
        _within(result.bodies["sphere"].velocity[:, axis],
                body[:, 7 + axis], speed_limit, f"{name} speed")
        # The source logs the reaction on the body with the opposite sign.
        # Its startup acceleration feedback is compared separately, so the
        # physical force gate begins after the first second.
        _within(result.body_forces["sphere"][100:, axis],
                -body[100:, 25 + axis], force_limit, f"{name} Morison load")
    for axis in (1, 3, 5):
        _within(result.bodies["sphere"].position[:, axis],
                body[:, 1 + axis], 1e-12, f"inactive position {axis}")
        _within(result.bodies["sphere"].velocity[:, axis],
                body[:, 7 + axis], 1e-12, f"inactive speed {axis}")
    assert np.max(np.abs(body[:, 1])) > 1
    assert np.max(np.abs(body[:, 5])) > .02

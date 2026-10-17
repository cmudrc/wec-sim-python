"""Pair a current-loaded normal/tangential Morison Sphere trajectory."""

import os
from pathlib import Path

import numpy as np

from wecsim import Current, RegularWave, WEC
from wecsim.hydrodynamics.morison import (
    MorisonElement, regular_morison_source_force,
    regular_wave_normal_tangential_morison_terms,
)


ROOT = Path(__file__).resolve().parents[1]
PORTABLE = ROOT / "tests/fixtures/matlab/sphere_normal_current/source.npz"
SOURCE_H5 = os.environ.get("WEC_SIM_SPHERE_H5")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
if bool(SOURCE_H5) != bool(REFERENCE):
    raise ValueError("live Morison current comparison needs both HDF5 and source output")
SPHERE_H5 = SOURCE_H5 or ROOT / "examples/data/sphere.h5"


def _source():
    """Prefer fresh MATLAB output, with the same full history pinned for CI."""
    if REFERENCE:
        folder = Path(REFERENCE)
        filenames = {
            "no_current_body": "SPHERE_MOVING_MORISON_NORMAL_1m-ME_body1.csv",
            "current_body": "SPHERE_MOVING_MORISON_NORMAL_CURRENT_1m-ME_body1.csv",
            "direct_delta": "SPHERE_MOVING_MORISON_NORMAL_CURRENT_direct_delta.csv",
            "current_wave": "SPHERE_MOVING_MORISON_NORMAL_CURRENT_wave.csv",
            "no_current_wave": "SPHERE_MOVING_MORISON_NORMAL_wave.csv",
        }
        return {key: np.loadtxt(folder / filename, delimiter=",")
                for key, filename in filenames.items()}
    with np.load(PORTABLE, allow_pickle=False) as saved:
        return {key: saved[key].copy() for key in saved.files}


def _response(current):
    wec = WEC("Sphere normal/tangential Morison current")
    sphere = wec.body(
        "sphere", Path(SPHERE_H5).resolve(),
        inertia=(20907301, 21306090.66, 37085481.11),
    )
    for axis in ("surge", "heave", "pitch"):
        wec.coordinate(axis, sphere.move(axis))
    wec.morison_element(
        sphere, point=sphere.at(0, 0, -2), axis=(0, 1, 0),
        drag_coefficient=(1, 0, 0),
        added_mass_coefficient=(1, 0, 0),
        area=(100, 0, 0), volume=20,
    )
    return wec.run(
        RegularWave(1, 8, current=current),
        dt=.01, end_time=40, ramp_time=10,
        initial_coordinate={"heave": 1},
    )


def _max_error(actual, expected, limit, label):
    error = float(np.max(np.abs(actual - expected)))
    assert error < limit, f"{label}: {error:.6g} exceeds {limit}"


def test_normal_current_source_force_and_independent_motion():
    source = _source()
    no_current = source["no_current_body"]
    current = source["current_body"]
    direct_delta = source["direct_delta"]
    wave = source["current_wave"]
    control_wave = source["no_current_wave"]
    assert no_current.shape == current.shape == (4001, 37)
    assert direct_delta.shape == (4001, 7)
    assert wave.shape == control_wave.shape == (4001, 2)
    assert all(np.isfinite(a).all() for a in source.values())
    np.testing.assert_allclose(wave[:, 0], np.arange(4001) * .01,
                               rtol=0, atol=1e-10)
    np.testing.assert_array_equal(wave, control_wave)

    element = MorisonElement(
        point=(0, 0, -2), axis=(0, 1, 0),
        drag_coefficient=(1, 0, 0),
        added_mass_coefficient=(1, 0, 0),
        area=(100, 0, 0), volume=20,
    )
    physical_increment = []
    source_increment = []
    for row in current:
        force_args = dict(
            time=row[0], position=row[1:7], velocity=row[7:13],
            wave_height=1, wave_period=8, ramp_time=10,
            water_depth=np.inf, rho=1000,
        )
        current_args = dict(
            current_speed=.8, current_direction=0,
            current_profile="power", current_depth=30,
        )
        physical_increment.append(
            regular_wave_normal_tangential_morison_terms(
                [element], **force_args, **current_args,
            )[0] - regular_wave_normal_tangential_morison_terms(
                [element], **force_args,
            )[0]
        )
        source_args = dict(
            **force_args, acceleration=row[31:37], direction=0,
            body_morison=2, element_axes=[(0, 1, 0)],
        )
        # Replace the physical infinite-depth label with the source's 200 m
        # HDF5 depth. Both use the same deep-water branch at this frequency.
        source_args["water_depth"] = 200
        source_increment.append(
            regular_morison_source_force(
                [element], **source_args, **current_args,
            ) - regular_morison_source_force([element], **source_args)
        )
    physical_increment = np.asarray(physical_increment)
    source_increment = np.asarray(source_increment)
    np.testing.assert_allclose(direct_delta[:, 0], current[:, 0], atol=1e-12)
    assert np.max(np.abs(source_increment - direct_delta[:, 1:])) < 2e-6
    for axis, minimum, label in ((0, 1000, "surge"),
                                 (2, 1000, "heave"),
                                 (4, 10000, "pitch")):
        peak = np.max(np.abs(direct_delta[:, 1 + axis]))
        assert peak > minimum, f"{label} current loading was not exercised"
        error = np.max(np.abs(
            physical_increment[:, axis] - direct_delta[:, 1 + axis]
        ))
        assert error < .02 * peak, (
            f"{label} physical/source current increment differs by {error:.6g}"
        )

    python_current = _response(Current(.8, 0, "power", 30))
    python_no_current = _response(None)
    _max_error(python_current.time, current[:, 0], 1e-10, "time")
    _max_error(python_current.wave_elevation, wave[:, 1], 1e-10,
               "wave elevation")
    for axis, label, absolute_limits, effect_minima, effect_limits in (
            (0, "surge", (.007, .0008), (5, .1), (.004, .001)),
            (2, "heave", (.002, .0025), (.01, .01), (.001, .001)),
            (4, "pitch", (.0003, .0002), (.005, .002), (.0003, .00015))):
        position_limit, speed_limit = absolute_limits
        effect_position_minimum, effect_speed_minimum = effect_minima
        effect_position_limit, effect_speed_limit = effect_limits
        actual = python_current.bodies["sphere"]
        baseline = python_no_current.bodies["sphere"]
        source_position_increment = (
            current[:, 1 + axis] - no_current[:, 1 + axis]
        )
        python_position_increment = (
            actual.position[:, axis] - baseline.position[:, axis]
        )
        source_speed_increment = (
            current[:, 7 + axis] - no_current[:, 7 + axis]
        )
        python_speed_increment = (
            actual.velocity[:, axis] - baseline.velocity[:, axis]
        )
        _max_error(actual.position[:, axis], current[:, 1 + axis],
                   position_limit, f"{label} position")
        _max_error(actual.velocity[:, axis], current[:, 7 + axis],
                   speed_limit, f"{label} speed")
        assert np.max(np.abs(source_position_increment)) > effect_position_minimum
        assert np.max(np.abs(source_speed_increment)) > effect_speed_minimum
        _max_error(python_position_increment, source_position_increment,
                   effect_position_limit, f"{label} current position response")
        _max_error(python_speed_increment, source_speed_increment,
                   effect_speed_limit, f"{label} current speed response")

"""Pair a wave-driven axial Morison element on the published Sphere model."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim import RegularWave, WEC
from wecsim.hydrodynamics.morison import (
    MorisonElement, regular_wave_axial_morison_terms,
    regular_wave_heave_morison_terms,
)


SPHERE_H5 = os.environ.get("WEC_SIM_SPHERE_H5")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
pytestmark = pytest.mark.skipif(
    not (SPHERE_H5 and REFERENCE),
    reason="paired wave-driven Morison Sphere output and HDF5 absent",
)


def _max_error(actual, expected, limit, label):
    assert actual.shape == expected.shape, label
    error = float(np.max(np.abs(actual - expected)))
    assert np.isfinite(error) and error < limit, (
        f"{label}: {error:.6g} exceeds {limit}"
    )


def test_wave_driven_sphere_morison_force_and_motion():
    folder = Path(REFERENCE)
    source = np.loadtxt(
        folder / "SPHERE_MOVING_MORISON_WAVE_1m-ME_body1.csv",
        delimiter=",",
    )
    wave = np.loadtxt(folder / "SPHERE_MOVING_MORISON_WAVE_wave.csv",
                      delimiter=",")
    source_force = np.loadtxt(
        folder / "SPHERE_MOVING_MORISON_WAVE_source_force.csv",
        delimiter=",",
    )
    assert source.shape == (4001, 37)
    assert wave.shape == (4001, 2)
    assert source_force.shape == (4001, 3)
    assert (np.isfinite(source).all() and np.isfinite(wave).all()
            and np.isfinite(source_force).all())
    element = MorisonElement(
        point=(0, 0, -2), drag_coefficient=(0, 0, 1),
        added_mass_coefficient=(0, 0, 1), area=(0, 0, 100), volume=20,
    )
    source_terms = np.array([
        regular_wave_heave_morison_terms(
            [element], center_z=-2, heave=row[3] + 2,
            speed=row[9], time=row[0], wave_height=1, wave_period=8,
            ramp_time=10, water_depth=np.inf, rho=1000,
        )
        for row in source
    ])
    source_physical = -source[:, 27]
    reconstructed = source_terms[:, 0] - source_terms[:, 1] * source[:, 33]
    _max_error(source_force[:, 0], source[:, 0], 1e-12,
               "source force replay time")
    # The source sphere's floating joint permits surge and pitch. Project
    # its saved state onto the supported heave coordinate for a strict
    # source-function comparison, then check the full Simulink force channel
    # against the replay made with its actual six-DOF state.
    _max_error(reconstructed, source_force[:, 2], 1e-6,
               "Morison force law at projected heave states")
    _max_error(source_force[50:, 1], source_physical[50:], 20,
               "full-state source force channel after startup")
    assert np.max(np.abs(source[:, 1])) > 1, "source surge was not exercised"

    wec = WEC("Sphere with wave-driven axial Morison element")
    sphere = wec.body("sphere", Path(SPHERE_H5).resolve())
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
    _max_error(result.time, source[:, 0], 1e-12, "time")
    _max_error(result.wave_elevation, wave[:, 1], 1e-10,
               "regular-wave elevation")
    _max_error(result.bodies["sphere"].position[:, 2], source[:, 3],
               .01, "sphere heave")
    _max_error(result.bodies["sphere"].velocity[:, 2], source[:, 9],
               .02, "sphere speed")
    physical = result.body_forces["sphere"][:, 2]
    _max_error(physical[100:], source_physical[100:], 1000,
               "wave-driven moving Morison force after startup")


def test_wave_driven_sphere_three_dof_morison_trajectory():
    source = np.loadtxt(
        Path(REFERENCE) / "SPHERE_MOVING_MORISON_WAVE_1m-ME_body1.csv",
        delimiter=",",
    )
    element = MorisonElement(
        point=(0, 0, -2), drag_coefficient=(0, 0, 1),
        added_mass_coefficient=(0, 0, 1), area=(0, 0, 100), volume=20,
    )
    source_force = -source[:, 25:31]
    # Check the physical proper-rotation force law on MATLAB's saved state
    # separately from the independently integrated trajectory. The pinned
    # source function uses a nonorthogonal general rotation matrix.
    on_source_state = np.stack([
        applied - added @ row[31:37]
        for row in source
        for applied, added in [regular_wave_axial_morison_terms(
            [element], position=row[1:7], velocity=row[7:13],
            time=row[0], wave_height=1, wave_period=8, ramp_time=10,
            water_depth=np.inf, rho=1000,
        )]
    ])
    for axis, limit, label in ((0, 500, "surge force"),
                               (2, 120, "heave force"),
                               (4, 220, "pitch moment")):
        _max_error(on_source_state[100:, axis], source_force[100:, axis],
                   limit, f"physical {label} on MATLAB state")

    wec = WEC("Sphere with three-DOF axial Morison element")
    sphere = wec.body(
        "sphere", Path(SPHERE_H5).resolve(),
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
    assert result.time.shape == (4001,)
    for axis, name, position_limit, velocity_limit in (
        (0, "surge", .015, .002),
        (2, "heave", .0015, .002),
        (4, "pitch", .00015, .00006),
    ):
        _max_error(result.bodies["sphere"].position[:, axis],
                   source[:, 1 + axis], position_limit, f"{name} position")
        _max_error(result.bodies["sphere"].velocity[:, axis],
                   source[:, 7 + axis], velocity_limit, f"{name} velocity")
    for axis, limit, label in ((0, 500, "surge force"),
                               (2, 200, "heave force"),
                               (4, 250, "pitch moment")):
        _max_error(result.body_forces["sphere"][100:, axis],
                   source_force[100:, axis], limit,
                   f"integrated {label} after startup")

"""Pair a derived moving normal/tangential Sphere case with MATLAB."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim import RegularWave, WEC
from wecsim.morison import (
    MorisonElement, regular_morison_source_force,
    regular_wave_normal_tangential_morison_terms,
)


SPHERE_H5 = os.environ.get("WEC_SIM_SPHERE_H5")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
pytestmark = pytest.mark.skipif(
    not (SPHERE_H5 and REFERENCE),
    reason="paired normal/tangential Sphere output and HDF5 absent",
)


def _max_error(actual, expected, limit, label):
    assert actual.shape == expected.shape, label
    assert np.isfinite(actual).all() and np.isfinite(expected).all(), label
    error = float(np.max(np.abs(actual - expected)))
    assert error < limit, f"{label}: {error:.6g} exceeds {limit}"


def test_sphere_moving_normal_tangential_trajectory():
    reference = Path(REFERENCE)
    source = np.loadtxt(
        reference / "SPHERE_MOVING_MORISON_NORMAL_1m-ME_body1.csv",
        delimiter=",",
    )
    source_wave = np.loadtxt(
        reference / "SPHERE_MOVING_MORISON_NORMAL_wave.csv",
        delimiter=",",
    )
    assert source.shape == (4001, 37)
    assert source_wave.shape == (4001, 2)

    wec = WEC("Sphere with moving normal/tangential Morison element")
    sphere = wec.body(
        "sphere", Path(SPHERE_H5).resolve(),
        inertia=(20907301, 21306090.66, 37085481.11),
    )
    for axis in ("surge", "heave", "pitch"):
        wec.coordinate(axis, sphere.move(axis))
    element = wec.morison_element(
        sphere, point=sphere.at(0, 0, -2), axis=(0, 1, 0),
        drag_coefficient=(1, 0, 0),
        added_mass_coefficient=(1, 0, 0),
        area=(100, 0, 0), volume=20,
    )
    result = wec.run(
        RegularWave(1, 8), dt=.01, end_time=40, ramp_time=10,
        initial_coordinate={"heave": 1},
    )
    assert result.case["bodies"][0]["morison_elements"][0]["axis"] == [0, 1, 0]
    _max_error(result.time, source[:, 0], 1e-10, "time")
    _max_error(result.wave_elevation, source_wave[:, 1], 1e-10,
               "wave elevation")
    for dof, name, position_limit, speed_limit, force_limit in (
        (0, "surge", .006, .0004, 50),
        (2, "heave", .0015, .0017, 150),
        (4, "pitch", .00025, .00012, 900),
    ):
        _max_error(result.bodies["sphere"].position[:, dof],
                   source[:, 1 + dof], position_limit, f"{name} position")
        _max_error(result.bodies["sphere"].velocity[:, dof],
                   source[:, 7 + dof], speed_limit, f"{name} speed")
        _max_error(result.body_forces["sphere"][100:, dof],
                   -source[100:, 25 + dof], force_limit,
                   f"{name} Morison force after startup")
    for dof in (1, 3, 5):
        _max_error(result.bodies["sphere"].position[:, dof],
                   source[:, 1 + dof], 1e-10,
                   f"inactive position DOF {dof}")
        _max_error(result.bodies["sphere"].velocity[:, dof],
                   source[:, 7 + dof], 1e-10,
                   f"inactive speed DOF {dof}")

    # The physical proper rotation and source option-2 rotation disagree at
    # nonzero pitch. Evaluate both on the same MATLAB states to separate
    # that force-law difference from independent trajectory drift.
    on_source_state = np.stack([
        applied - mass @ row[31:37]
        for row in source
        for applied, mass in [regular_wave_normal_tangential_morison_terms(
            [element], position=row[1:7], velocity=row[7:13],
            time=row[0], wave_height=1, wave_period=8,
            water_depth=np.inf, ramp_time=10, rho=1000,
        )]
    ])
    for dof, limit in ((0, 50), (2, 80), (4, 900)):
        _max_error(on_source_state[100:, dof],
                   -source[100:, 25 + dof], limit,
                   f"physical force on source states DOF {dof}")

    # The source-convention diagnostic reconstructs the logged force much
    # more closely; the remaining physical gap comes from its different
    # rotation and point-velocity kinematics. At this 8 s period, 200 m
    # uses the same deep-water branch as the source's infinite depth.
    source_replay = np.stack([
        regular_morison_source_force(
            [element], time=row[0], position=row[1:7],
            velocity=row[7:13], acceleration=row[31:37],
            wave_height=1, wave_period=8, direction=0,
            water_depth=200, ramp_time=10, rho=1000,
            body_morison=2, element_axes=[(0, 1, 0)],
        ) for row in source
    ])
    for dof, limit in ((0, 2), (2, 12), (4, 3)):
        _max_error(source_replay[100:, dof],
                   -source[100:, 25 + dof], limit,
                   f"source-convention force replay DOF {dof}")

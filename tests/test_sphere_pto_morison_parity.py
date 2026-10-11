"""Pair one moving Sphere Morison element and a linear PTO with MATLAB."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim import RegularWave, WEC, WorldPoint


SPHERE_H5 = os.environ.get("WEC_SIM_SPHERE_H5")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
pytestmark = pytest.mark.skipif(
    not (SPHERE_H5 and REFERENCE),
    reason="paired Sphere PTO and Morison MATLAB output is absent",
)


def _check(actual, expected, limit, label):
    actual = np.asarray(actual)
    expected = np.asarray(expected)
    assert actual.shape == expected.shape, label
    assert np.isfinite(actual).all() and np.isfinite(expected).all(), label
    error = float(np.max(np.abs(actual - expected)))
    assert error < limit, f"{label}: {error:.6g} exceeds {limit}"


def test_moving_normal_tangential_morison_with_linear_pto():
    source = Path(REFERENCE)
    body = np.loadtxt(
        source / "SPHERE_PTO_MORISON_NORMAL_ConfiguredMorison_body1.csv",
        delimiter=",",
    )
    pto = np.loadtxt(
        source / "SPHERE_PTO_MORISON_NORMAL_ConfiguredMorison_pto1.csv",
        delimiter=",",
    )
    controller = np.loadtxt(
        source / "SPHERE_PTO_MORISON_NORMAL_controller.csv", delimiter=",",
    )
    wave = np.loadtxt(
        source / "SPHERE_PTO_MORISON_NORMAL_wave.csv", delimiter=",",
    )
    assert body.shape == (2001, 37)
    assert pto.shape == (2001, 25)
    assert controller.shape == (2001, 13)
    assert wave.shape == (2001, 2)

    wec = WEC("Sphere PTO and moving normal/tangential Morison")
    sphere = wec.body(
        "sphere", Path(SPHERE_H5).resolve(),
        inertia=(20907301, 21306090.66, 37085481.11),
    )
    wec.coordinate("heave", sphere.move("heave"))
    wec.morison_element(
        sphere, point=sphere.at(0, 0, -2), axis=(0, 1, 0),
        drag_coefficient=(1, 0, 0), added_mass_coefficient=(1, 0, 0),
        area=(100, 0, 0), volume=20,
    )
    wec.pto(
        "main", sphere.at(1, 0, 2), WorldPoint(1, 0, 0),
        axis=(0, 0, 1), stiffness=50_000, damping=960_870,
    )
    result = wec.run(
        RegularWave(2.5, 9.6664), dt=.02, end_time=40, ramp_time=10,
    )
    _check(result.time, body[:, 0], 1e-10, "body time")
    _check(result.time, pto[:, 0], 1e-10, "PTO time")
    _check(result.time, controller[:, 0], 1e-10, "controller time")
    _check(result.wave_elevation, wave[:, 1], 1e-10, "wave elevation")
    _check(result.bodies["sphere"].position[:, 2], body[:, 3], .005,
           "heave position")
    _check(result.bodies["sphere"].velocity[:, 2], body[:, 9], .005,
           "heave speed")
    _check(-result.ptos["main"].stroke, pto[:, 3], .005, "PTO stroke")
    _check(-result.ptos["main"].velocity, pto[:, 9], .005, "PTO speed")
    _check(-result.ptos["main"].force, pto[:, 15] + controller[:, 3], 5000,
           "combined PTO force")
    _check(result.body_forces["sphere"][100:, 2], -body[100:, 27], 5000,
           "physical heave Morison force after startup")

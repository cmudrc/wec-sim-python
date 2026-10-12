"""Pair one moving Sphere Morison element and a linear PTO with MATLAB."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim import RegularWave, WEC, WorldPoint
from wecsim.hydrodynamics.morison import (
    MorisonElement, regular_wave_normal_tangential_morison_terms,
)


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
        axis=(0, 0, 1), stiffness=50_000, damping=100_000,
        velocity_control_gain=860_870,
    )
    result = wec.run(
        RegularWave(2.5, 9.6664), dt=.02, end_time=40, ramp_time=10,
    )
    _check(result.time, body[:, 0], 1e-10, "body time")
    _check(result.time, pto[:, 0], 1e-10, "PTO time")
    _check(result.time, controller[:, 0], 1e-10, "controller time")
    _check(result.time, wave[:, 0], 1e-10, "wave time")
    _check(result.wave_elevation, wave[:, 1], 1e-10, "wave elevation")
    for dof in (0, 1, 3, 4, 5):
        _check(result.bodies["sphere"].position[:, dof], body[:, 1 + dof],
               1e-10, f"inactive position DOF {dof}")
        _check(result.bodies["sphere"].velocity[:, dof], body[:, 7 + dof],
               1e-10, f"inactive speed DOF {dof}")
    _check(result.bodies["sphere"].position[:, 2], body[:, 3], .0002,
           "heave position")
    _check(result.bodies["sphere"].velocity[:, 2], body[:, 9], .0002,
           "heave speed")
    _check(-result.ptos["main"].stroke, pto[:, 3], .0002, "PTO stroke")
    _check(-result.ptos["main"].velocity, pto[:, 9], .0002, "PTO speed")
    assert result.case["ptos"][0]["damping"] == 100_000
    assert result.case["ptos"][0]["velocity_control_gain"] == 860_870
    native_force = -(50_000 * result.ptos["main"].stroke
                     + 100_000 * result.ptos["main"].velocity)
    controller_force = result.ptos["main"].control_force
    controller_power = result.ptos["main"].control_absorbed_power
    assert controller_force is not None and controller_power is not None
    np.testing.assert_allclose(
        controller_force, -860_870 * result.ptos["main"].velocity,
        rtol=0, atol=1e-10,
    )
    np.testing.assert_allclose(
        result.ptos["main"].force, native_force + controller_force,
        rtol=0, atol=1e-9,
    )
    _check(-native_force, pto[:, 15], 20, "native PTO force")
    _check(-controller_force, controller[:, 3], 100,
           "proportional controller force")
    _check(-result.ptos["main"].force, pto[:, 15] + controller[:, 3], 120,
           "combined PTO force")
    _check(native_force * result.ptos["main"].velocity, pto[:, 21], 10,
           "native PTO power")
    _check(-controller_power, controller[:, 9], 50, "controller power")
    for dof, name, limit in (
        (0, "surge", 5), (2, "heave", 20), (4, "pitch", 10),
    ):
        _check(result.body_forces["sphere"][:, dof], -body[:, 25 + dof],
               limit, f"{name} Morison force")

    # On the saved MATLAB motion, the proper physical option-2 law matches
    # the source force directly because this PTO joint holds the body upright.
    element = MorisonElement(
        (0, 0, -2), (1, 0, 0), (1, 0, 0), (100, 0, 0), 20,
        axis=(0, 1, 0),
    )
    on_source_state = np.stack([
        applied - mass @ row[31:37]
        for row in body
        for applied, mass in [regular_wave_normal_tangential_morison_terms(
            [element], position=row[1:7], velocity=row[7:13],
            time=row[0], wave_height=2.5, wave_period=9.6664,
            water_depth=np.inf, ramp_time=10, rho=1000,
        )]
    ])
    _check(on_source_state[:, 0], -body[:, 25], 1e-6,
           "surge Morison law on source motion")
    _check(on_source_state[:, 2], -body[:, 27], 5,
           "heave Morison law on source motion")
    _check(on_source_state[:, 4], -body[:, 29], 1e-6,
           "pitch Morison law on source motion")

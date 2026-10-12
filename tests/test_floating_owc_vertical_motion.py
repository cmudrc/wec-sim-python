"""Pair independent two-body heave motion under saved external OWC loads."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim.hydrodynamics.bodyClass import BodyClass
from wecsim.dynamics.generalDynamics import BodyMotion, DynamicBody, GeneralizedDynamics


REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
HYDRO = os.environ.get("WEC_SIM_FLOATING_OWC_H5")
pytestmark = pytest.mark.skipif(
    not (REFERENCE and HYDRO),
    reason="pinned floating OWC bodies, loads, and HDF5 absent",
)


def _read(name):
    return np.loadtxt(Path(REFERENCE) / f"FLOATING_OWC_SOURCE_{name}.csv",
                      delimiter=",", ndmin=2)


def _body(number, mooring, pressure):
    body = BodyClass(HYDRO)
    body.bodyNumber = number
    body.bodyTotal = 2
    body.readH5file()
    column_mass = 4493450.0
    if number == 1:
        body.mass = "equilibrium"
        inertia = np.array([1.531e9, 1.531e9, 0.1118e9])
    else:
        body.mass = column_mass
        radius = 5.89 / 2
        transverse = column_mass * (3 * radius**2 + 50.69**2) / 12
        inertia = np.array([transverse, transverse,
                            0.5 * column_mass * radius**2])
    body.inertia = inertia
    body.hydroStiffness = np.zeros((6, 6))
    body.viscDrag = {
        "Drag": np.zeros((6, 6)), "cd": np.zeros(6),
        "characteristicArea": np.zeros(6),
    }
    body.linearDamping = np.zeros((6, 6))
    frequency = 2 * np.pi / 11.2
    body.hydroForcePre(
        frequency, [0], 1, np.array([0.0]), [], 0.01, 1000, 9.81,
        "regular", np.zeros((2, 2)), number, 2, 0, 0, 0,
    )
    mass = float(np.asarray(body.mass).item())
    center = np.asarray(body.cg, dtype=float).ravel()
    mapping = np.zeros((6, 2))
    mapping[2, number - 1] = 1

    def motion(coordinate, speed):
        return BodyMotion(mapping @ coordinate, mapping, np.zeros(6))

    real = np.asarray(body.hydroForce["fExt"]["re"])
    imaginary = np.asarray(body.hydroForce["fExt"]["im"])
    area = np.pi * 5.89**2 / 4

    def excitation(time, coordinate, speed):
        ramp = 1 if time >= 50 else (1 - np.cos(np.pi * time / 50)) / 2
        force = 2.25 * ramp * (
            real * np.cos(frequency * time)
            - imaginary * np.sin(frequency * time)
        )
        chamber_force = area * np.interp(time, pressure[:, 0], pressure[:, 1])
        if number == 1:
            force[2] += (np.interp(time, mooring[:, 0], mooring[:, 15])
                         + chamber_force)
        else:
            force[2] -= chamber_force
        return force

    added = [np.zeros((6, 6)) for _ in range(2)]
    damping = [np.zeros((6, 6)) for _ in range(2)]
    added[number - 1] = np.asarray(body.hydroForce["fAddedMass"])
    damping[number - 1] = np.asarray(body.hydroForce["fDamping"])
    return DynamicBody(
        rigid_mass=np.diag([mass] * 3 + inertia.tolist()),
        added_mass=tuple(added), damping=tuple(damping),
        restoring=np.asarray(body.hydroForce["linearHydroRestCoef"]),
        static_force=np.array([
            0, 0, (1000 * float(np.asarray(body.dispVol).item()) - mass)
            * 9.81, 0, 0, 0,
        ]),
        reference_position=np.r_[center, np.zeros(3)],
        motion=motion, excitation=lambda time: np.zeros(6),
        state_excitation=excitation,
    )


def test_floating_owc_heave_dynamics_with_prescribed_mooring_and_pressure():
    source = [_read(f"FloatingOWC_body{number}") for number in (1, 2)]
    mooring = _read("coupling")
    pressure = _read("x_deltaP")
    assert len(source[0]) == len(source[1]) == len(mooring) == len(pressure) == 50001
    assert np.min(mooring[:, 15]) < -1e5
    assert np.max(np.abs(pressure[:, 1])) > 1000

    system = GeneralizedDynamics(tuple(
        _body(number, mooring, pressure) for number in (1, 2)
    ), 2)
    response = system.integrate(dt=0.1, end_time=500,
                                adaptive_regular=True)
    for index, body in enumerate(source):
        sampled = body[::10]
        np.testing.assert_allclose(response.time, sampled[:, 0],
                                   rtol=0, atol=1e-9)
        heave_error = response.coordinate[:, index] - (
            sampled[:, 3] - body[0, 3]
        )
        speed_error = response.speed[:, index] - sampled[:, 9]
        assert np.max(np.abs(heave_error)) < 0.015
        assert np.max(np.abs(speed_error)) < 0.006

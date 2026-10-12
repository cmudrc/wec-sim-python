"""Check numeric output from the pinned floating OWC MATLAB application."""

import os
from pathlib import Path

import h5py
import numpy as np
import pytest

from wecsim.hydrodynamics.bodyClass import BodyClass


REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
HYDRO = os.environ.get("WEC_SIM_FLOATING_OWC_H5")
pytestmark = pytest.mark.skipif(
    not (REFERENCE and HYDRO),
    reason="pinned floating OWC MATLAB output and BEMIO HDF5 not provided",
)


def _read(name):
    values = np.loadtxt(Path(REFERENCE) / f"FLOATING_OWC_SOURCE_{name}.csv",
                        delimiter=",")
    assert values.ndim == 2 and np.isfinite(values).all()
    assert np.all(np.diff(values[:, 0]) > 0)
    return values


def _static(name):
    values = np.loadtxt(Path(REFERENCE) / f"FLOATING_OWC_SOURCE_{name}.csv",
                        delimiter=",", ndmin=2)
    assert np.isfinite(values).all()
    return values


def _prepared_body(number):
    body = BodyClass(HYDRO)
    body.bodyNumber = number
    body.bodyTotal = 2
    body.readH5file()
    body.mass = "equilibrium" if number == 1 else 4493450
    body.inertia = np.zeros(3)
    body.hydroStiffness = np.zeros((6, 6))
    body.viscDrag = {
        "Drag": np.zeros((6, 6)), "cd": np.zeros(6),
        "characteristicArea": np.zeros(6),
    }
    body.linearDamping = np.zeros((6, 6))
    body.hydroForcePre(
        2 * np.pi / 11.2, [0], 1, np.array([0.0]), [], 0.01,
        1000, 9.81, "regular", np.zeros((2, 2)), number, 2, 0, 0, 0,
    )
    return body


def test_published_floating_owc_source_has_coupled_motion_and_power():
    wave = _read("wave")
    assert wave.shape[1] == 2
    assert abs(wave[-1, 0] - 500) < 1e-8
    assert np.max(np.abs(wave[:, 1])) > 1

    bodies = [_read(f"FloatingOWC_body{number}") for number in (1, 2)]
    for body in bodies:
        assert body.shape[1] == 25
        np.testing.assert_allclose(body[:, 0], wave[:, 0], rtol=0, atol=1e-8)
    assert np.ptp(bodies[0][:, 3] - bodies[1][:, 3]) > 0.01

    coupling = _read("coupling")
    assert coupling.shape[1] == 19
    assert np.max(np.linalg.norm(coupling[:, 13:16], axis=1)) > 100
    fairlead = _read("fairlead_tension")
    assert fairlead.shape[1] >= 4
    assert np.max(fairlead[:, 1:]) > 100

    for name in ("x_xOWC", "x_signal4", "x_deltaP", "x_vTurb", "u",
                 "P_pneumatic", "P_turb", "eta_turb"):
        trace = _read(name)
        assert trace.shape[1] >= 2
        assert trace[-1, 0] >= 499.9
    with h5py.File(HYDRO) as hydro:
        assert "body1" in hydro and "body2" in hydro
        assert "hydro_coeffs" in hydro["body1"]


def test_floating_owc_python_regular_excitation_matches_both_matlab_bodies():
    """Pair the production HDF5 wave-force path with the published source run."""
    period = 11.2
    frequency = 2 * np.pi / period
    wave = _read("wave")
    time = wave[:, 0]
    ramp = np.ones_like(time)
    early = time < 50
    ramp[early] = (1 - np.cos(np.pi * time[early] / 50)) / 2
    amplitude = 4.5 / 2
    np.testing.assert_allclose(
        wave[:, 1], amplitude * ramp * np.cos(frequency * time),
        rtol=1e-11, atol=1e-11,
    )

    for number in (1, 2):
        body = _prepared_body(number)
        coefficients = body.hydroForce["fExt"]
        predicted = amplitude * ramp[:, None] * (
            np.cos(frequency * time)[:, None] * coefficients["re"]
            - np.sin(frequency * time)[:, None] * coefficients["im"]
        )
        source = _read(f"FloatingOWC_body{number}")
        np.testing.assert_allclose(source[:, 0], time, rtol=0, atol=1e-8)
        excitation = source[:, 19:25]
        assert np.max(np.abs(excitation)) > 1e5
        assert np.max(np.abs(predicted - excitation)) < 1e-5


def test_floating_owc_source_exports_mechanical_force_balance():
    for number in (1, 2):
        source = _read(f"FloatingOWC_body{number}")
        forces = _read(f"body{number}_forces")
        assert forces.shape == (5001, 37)
        np.testing.assert_allclose(forces[:, 0], source[::10, 0],
                                   rtol=0, atol=1e-8)
        mass = _static(f"body{number}_mass")
        assert mass.shape == (1, 4) and mass[0, 0] > 0
        # The published input leaves b2b unset: each body uses its own 6x6
        # hydrodynamic blocks despite the shared 12-DOF BEMIO file.
        for name, shape in (("added_mass", (6, 6)),
                            ("radiation_damping", (6, 6)),
                            ("hydrostatic", (6, 6))):
            matrix = _static(f"body{number}_{name}")
            assert matrix.shape == shape
        body = _prepared_body(number)
        np.testing.assert_allclose(
            _static(f"body{number}_added_mass"),
            body.hydroForce["fAddedMass"], rtol=1e-12, atol=1e-5,
        )
        damping = _static(f"body{number}_radiation_damping")
        np.testing.assert_allclose(
            damping, body.hydroForce["fDamping"],
            rtol=1e-12, atol=1e-7,
        )
        np.testing.assert_allclose(
            _static(f"body{number}_hydrostatic"),
            body.hydroForce["linearHydroRestCoef"],
            rtol=1e-12, atol=1e-5,
        )
        np.testing.assert_allclose(
            forces[:, 1:7], source[::10, 7:13] @ damping.T,
            rtol=1e-12, atol=1e-8,
        )

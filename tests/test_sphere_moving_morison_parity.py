"""Paired Sphere free decay with an active heave Morison element."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim import NoWave, WEC
from wecsim.hydrodynamics.morison import MorisonElement, no_wave_heave_morison_terms


SPHERE_H5 = os.environ.get("WEC_SIM_SPHERE_H5")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
pytestmark = pytest.mark.skipif(
    not (SPHERE_H5 and REFERENCE),
    reason="paired Sphere moving-Morison output and HDF5 absent",
)


def _max_error(actual, expected, limit, label):
    assert actual.shape == expected.shape, label
    error = float(np.max(np.abs(actual - expected)))
    assert np.isfinite(error) and error < limit, (
        f"{label}: {error:.6g} exceeds {limit}"
    )


def test_sphere_heave_morison_force_and_trajectory():
    source = np.loadtxt(
        Path(REFERENCE) / "SPHERE_MOVING_MORISON_1m-ME_body1.csv",
        delimiter=",",
    )
    assert source.shape == (4001, 37)
    assert np.isfinite(source).all()
    element = MorisonElement(
        point=(0, 0, -2), drag_coefficient=(0, 0, 1),
        added_mass_coefficient=(0, 0, 1), area=(0, 0, 100), volume=20,
    )
    # The source output stores forceMorisonAndViscous with the opposite sign
    # to the physical force on the body, as in the fixed Morison case.
    expected_physical = -source[:, 27]
    terms = np.array([
        no_wave_heave_morison_terms(
            [element], center_z=-2, heave=position + 2,
            speed=speed, rho=1000,
        )
        for position, speed in zip(source[:, 3], source[:, 9])
    ])
    reconstructed = terms[:, 0] - terms[:, 1] * source[:, 33]
    # MATLAB's acceleration feedback is not the saved output acceleration
    # during startup. Once settled, the saved acceleration reconstructs the
    # applied source force closely; the early mismatch remains diagnostic.
    _max_error(reconstructed[50:], expected_physical[50:], 20,
               "Morison force on MATLAB's saved state after startup")
    assert np.max(np.abs(reconstructed[:50] - expected_physical[:50])) > 30_000
    assert np.max(np.abs(expected_physical)) > 1000

    wec = WEC("Sphere with active heave Morison element")
    sphere = wec.body("sphere", Path(SPHERE_H5).resolve())
    wec.coordinate("heave", sphere.move("heave"))
    wec.morison_element(
        sphere, point=sphere.at(0, 0, -2),
        drag_coefficient=(0, 0, 1),
        added_mass_coefficient=(0, 0, 1),
        area=(0, 0, 100), volume=20,
    )
    result = wec.run(
        NoWave(), dt=.01, end_time=40, radiation_memory=15,
        initial_coordinate={"heave": 1},
    )
    _max_error(result.time, source[:, 0], 1e-12, "time")
    _max_error(result.bodies["sphere"].position[:, 2], source[:, 3],
               .0005, "sphere heave")
    _max_error(result.bodies["sphere"].velocity[:, 2], source[:, 9],
               .0011, "sphere speed")
    physical_force = result.body_forces["sphere"][:, 2]
    _max_error(physical_force[100:], expected_physical[100:],
               50, "moving Morison heave force after startup")
    impulse_difference = abs(np.trapezoid(
        physical_force - expected_physical, result.time,
    ))
    assert impulse_difference < 250, (
        f"moving Morison force impulse differs by {impulse_difference:.6g} N s"
    )

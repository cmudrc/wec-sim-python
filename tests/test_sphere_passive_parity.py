"""Compare the published Sphere passive-controller case with Python dynamics."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim import RegularWave, WEC, WorldPoint


SPHERE_H5 = os.environ.get("WEC_SIM_SPHERE_H5")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
ROOT = Path(__file__).resolve().parents[1]
PORTABLE = Path(__file__).parent / "fixtures/matlab/sphere_passive.npz"


def _source():
    if SPHERE_H5 and REFERENCE:
        reference = Path(REFERENCE)
        body = np.loadtxt(reference / "Sphere_Passive_Passive_P_body1.csv", delimiter=",")
        controller = np.loadtxt(reference / "Sphere_Passive_controller.csv", delimiter=",")
        assert body.shape == (20_001, 25)
        assert controller.shape == (20_001, 13)
        np.testing.assert_allclose(body[:, 0], controller[:, 0], rtol=0, atol=1e-10)
        return (body[:, 0], body[:, 3], body[:, 9],
                controller[:, 3], controller[:, 9], Path(SPHERE_H5).resolve())
    if not SPHERE_H5 and not REFERENCE:
        with np.load(PORTABLE) as saved:
            channels = tuple(saved[key] for key in
                             ("time", "position", "velocity", "force", "power"))
        return (*channels, ROOT / "examples/data/sphere.h5")
    pytest.fail("provide both WEC_SIM_SPHERE_H5 and WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")


def test_sphere_passive_controller_against_matlab():
    time, position, velocity, source_force, source_power, hydro = _source()
    assert time.shape == (20_001,)
    wec = WEC("Sphere passive controller")
    sphere = wec.body(
        "sphere", hydro,
        inertia=(20_907_301, 21_306_090.66, 37_085_481.11),
    )
    wec.coordinate("heave", sphere.move("heave"))
    wec.pto("passive", sphere.at(0, 0, 2), WorldPoint(0, 0, 0),
            axis=(0, 0, 1), damping=860_870)
    response = wec.run(
        RegularWave(2.5, 9.6664), dt=0.02, end_time=400, ramp_time=10,
    )
    np.testing.assert_allclose(response.time, time, rtol=0, atol=1e-10)
    python_position = response.bodies["sphere"].position[:, 2]
    python_velocity = response.bodies["sphere"].velocity[:, 2]
    # The Python PTO reports force along its from-to axis and positive
    # absorbed power. MATLAB's controller channels use the opposite signs.
    force = -response.ptos["passive"].force
    power = -response.ptos["passive"].absorbed_power
    for name, actual, expected, limit in (
        ("position", python_position, position, 1.5e-4),
        ("velocity", python_velocity, velocity, 1.5e-4),
        ("controller force", force, source_force, 150),
        ("controller power", power, source_power, 100),
    ):
        assert actual.shape == expected.shape, name
        assert np.isfinite(expected).all(), name
        error = np.max(np.abs(actual - expected))
        assert error < limit, f"{name}: max error {error:.6g} exceeds {limit}"

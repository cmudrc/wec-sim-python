"""Pair core geometry force preprocessing with pinned MATLAB bodyClass."""

import os
from pathlib import Path

import numpy as np
import pytest
from scipy.io import loadmat

from wecsim.hydrodynamics.bodyClass import BodyClass


ROOT = Path(__file__).resolve().parents[1]
FRESH = os.environ.get("WEC_SIM_MATLAB_REFERENCE_DIR")
CORE = os.environ.get("WEC_SIM_MATLAB_CORE_DIR")
FIXTURE = ROOT / "tests/fixtures/matlab/core_hydro_preprocess.npz"
FIELDS = ("linearRestoring", "addedMass", "radiationKernel",
          "excitationReal", "excitationImag")


def _source(geometry, number):
    if FRESH:
        values = loadmat(Path(FRESH) / f"core_hydro_{geometry}_body{number}.mat")
        return {field: values[field] for field in FIELDS}
    with np.load(FIXTURE) as values:
        return {field: values[f"{geometry}_{number}_{field}"]
                for field in FIELDS}


def _configured_body(geometry, number):
    if geometry == "rm3_current":
        file = Path(CORE) / "examples/RM3/hydroData/rm3.h5"
    else:
        file = ROOT / "examples/data" / f"{geometry}.h5"
    body = BodyClass(str(file))
    body.bodyNumber = number
    body.bodyTotal = 2
    body.readH5file()
    body.hydroStiffness = np.zeros((6, 6))
    body.viscDrag = {"Drag": np.zeros((6, 6)),
                     "cd": np.zeros(6), "characteristicArea": np.zeros(6)}
    body.linearDamping = np.zeros((6, 6))
    body.mass = "equilibrium"
    return body


@pytest.mark.parametrize("geometry", ("rm3", "oswec", "rm3_current"))
@pytest.mark.parametrize("number", (1, 2))
def test_regular_cic_force_preprocessing_against_matlab(geometry, number):
    """Check hydrostatics, radiation memory and excitation on both bodies."""
    if geometry == "rm3_current" and not (FRESH and CORE):
        pytest.skip("current core RM3 HDF5 and fresh MATLAB source needed")
    body = _configured_body(geometry, number)
    time = np.arange(101) * 0.1
    body.hydroForcePre(
        2 * np.pi / 8, [0], len(time), time, [], 0.1, 1000, 9.81,
        "regularCIC", np.vstack((time, np.zeros_like(time))),
        number, [], 0, 0, 0,
    )
    force = body.hydroForce
    expected = _source(geometry, number)
    actual = {
        "linearRestoring": force["linearHydroRestCoef"],
        "addedMass": force["fAddedMass"],
        "radiationKernel": force["irkb"],
        "excitationReal": force["fExt"]["re"],
        "excitationImag": force["fExt"]["im"],
    }
    for field in FIELDS:
        source = np.asarray(expected[field])
        python = np.asarray(actual[field])
        if field.startswith("excitation"):
            assert source.shape == (1, 6)
            assert python.shape == (6,)
            source = source[0]
        else:
            assert python.shape == source.shape, (geometry, number, field)
        assert np.isfinite(source).all() and np.max(np.abs(source)) > 0
        np.testing.assert_allclose(
            python, source, rtol=1e-10, atol=1e-5,
            err_msg=f"{geometry} body {number} {field}",
        )


@pytest.mark.parametrize("number", (1, 2))
def test_oswec_irregular_directional_excitation_against_matlab(number):
    """Check 64 frequency bins and three incident headings on each body."""
    body = _configured_body("oswec", number)
    bem_frequency = np.asarray(body.hydroData["simulation_parameters"]["w"]).ravel()
    omega = np.linspace(bem_frequency.min(), bem_frequency.max(), 64)
    time = np.arange(101) * 0.1
    body.hydroForcePre(
        omega, [0, 30, 90], len(time), time, 64, 0.1, 1000, 9.81,
        "irregular", np.vstack((time, np.zeros_like(time))),
        number, [], 0, 0, 0,
    )
    if FRESH:
        expected = loadmat(Path(FRESH) / f"core_hydro_oswec_irregular_body{number}.mat")
    else:
        with np.load(FIXTURE) as values:
            expected = {field: values[f"oswec_irregular_{number}_{field}"]
                        for field in ("omega", "excitationReal", "excitationImag")}
    np.testing.assert_allclose(omega, expected["omega"].ravel(), rtol=0, atol=1e-12)
    for field, key in (("excitationReal", "re"), ("excitationImag", "im")):
        source = np.asarray(expected[field])
        python = np.asarray(body.hydroForce["fExt"][key])
        assert source.shape == python.shape == (3, 64, 6)
        assert np.isfinite(source).all() and np.max(np.abs(source)) > 0
        np.testing.assert_allclose(
            python, source, rtol=1e-10, atol=1e-5,
            err_msg=f"OSWEC body {number} {field}",
        )

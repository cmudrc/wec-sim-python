"""Pair core geometry force preprocessing with pinned MATLAB bodyClass."""

from copy import deepcopy
import hashlib
import os
from pathlib import Path

import numpy as np
import pytest
from scipy.io import loadmat, savemat

from wecsim import ImportedElevationWave, WEC
from wecsim.hydrodynamics.bodyClass import BodyClass


ROOT = Path(__file__).resolve().parents[1]
FRESH = os.environ.get("WEC_SIM_MATLAB_REFERENCE_DIR")
CORE = os.environ.get("WEC_SIM_MATLAB_CORE_DIR")
CURRENT_RM3_SHA256 = (
    "a909360c6479773c113f792be72afb13"
    "51cae9299d62e501fa80561dca450aa3"
)
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


def _configured_body(geometry, number, source="hdf5"):
    if geometry == "rm3_current":
        file = Path(CORE) / "examples/RM3/hydroData/rm3.h5"
        assert hashlib.sha256(file.read_bytes()).hexdigest() == CURRENT_RM3_SHA256
    else:
        file = ROOT / "examples/data" / f"{geometry}.h5"
    body = BodyClass(str(file))
    body.bodyNumber = number
    body.bodyTotal = 2
    body.readH5file()
    if source == "memory":
        # MATLAB's MCR route supplies one body's BEM structure directly.
        loaded = BodyClass(None)
        loaded.bodyNumber = number
        loaded.bodyTotal = 2
        loaded.loadHydroData(deepcopy(body.hydroData))
        body = loaded
    body.hydroStiffness = np.zeros((6, 6))
    body.viscDrag = {"Drag": np.zeros((6, 6)),
                     "cd": np.zeros(6), "characteristicArea": np.zeros(6)}
    body.linearDamping = np.zeros((6, 6))
    body.mass = "equilibrium"
    return body


@pytest.mark.parametrize("geometry", ("rm3", "oswec", "rm3_current"))
@pytest.mark.parametrize("number", (1, 2))
@pytest.mark.parametrize("source", ("hdf5", "memory"))
def test_regular_cic_force_preprocessing_against_matlab(geometry, number, source):
    """Check hydrostatics, radiation memory and excitation on both bodies."""
    if geometry == "rm3_current" and not (FRESH and CORE):
        pytest.skip("current core RM3 HDF5 and fresh MATLAB source needed")
    body = _configured_body(geometry, number, source)
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
@pytest.mark.parametrize("source", ("hdf5", "memory"))
def test_oswec_irregular_directional_excitation_against_matlab(number, source):
    """Check 64 frequency bins and three incident headings on each body."""
    body = _configured_body("oswec", number, source)
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


def _imported_heading_source(number):
    if FRESH:
        values = loadmat(Path(FRESH) / f"core_hydro_oswec_heading_body{number}.mat")
        return values["waveAmpTime"], values["excitationForce"]
    with np.load(FIXTURE) as values:
        return (values[f"oswec_heading_{number}_waveAmpTime"],
                values[f"oswec_heading_{number}_excitationForce"])


@pytest.mark.parametrize("number", (1, 2))
@pytest.mark.parametrize("source", ("hdf5", "memory"))
def test_oswec_imported_elevation_heading_against_matlab(number, source, tmp_path):
    """An interior BEM heading changes all six imported-wave force channels."""
    time = np.arange(201) * .1
    eta = .4 * np.cos(2 * np.pi * time / 7) + .15 * np.sin(2 * np.pi * time / 3)
    record = tmp_path / "etaData.mat"
    savemat(record, {"etaData": np.column_stack((time, eta))})
    source_wave, source_force = _imported_heading_source(number)
    assert source_wave.shape == (201, 2)
    assert source_force.shape == (201, 6)
    ramp = np.ones_like(time)
    early = time < 2
    ramp[early] = (1 - np.cos(np.pi * time[early] / 2)) / 2
    expected_wave = eta * ramp
    np.testing.assert_allclose(source_wave[:, 0], time, rtol=0, atol=1e-12)
    np.testing.assert_allclose(source_wave[:, 1], expected_wave, rtol=0, atol=1e-12)

    body = _configured_body("oswec", number, source)
    body.hydroForce["userDefinedFe"] = np.zeros((len(time), 6))
    body.userDefinedExcitation(
        np.vstack((time, expected_wave)), .1, [15], 1000, 9.81,
    )
    np.testing.assert_allclose(body.hydroForce["userDefinedFe"], source_force,
                               rtol=1e-10, atol=1e-5)
    with pytest.raises(ValueError, match="outside the hydro database headings"):
        body.userDefinedExcitation(
            np.vstack((time, expected_wave)), .1, [95], 1000, 9.81,
        )

    if number == 1 and source == "hdf5":
        # The public wave setting must reach production body preprocessing,
        # not just the direct BodyClass call above.
        wec = WEC("OSWEC imported heading")
        flap = wec.body("flap", ROOT / "examples/data/oswec.h5",
                        inertia=(1e7, 1e7, 1e7))
        wec.coordinate("heave", flap.move("heave"))
        result = wec.run(
            ImportedElevationWave(record, direction=15),
            dt=.1, end_time=20, ramp_time=2, radiation_memory=10,
        )
        np.testing.assert_allclose(result.wave_elevation, expected_wave,
                                   rtol=0, atol=1e-12)
        force = dict(result.raw.extra_outputs)["body1_excitation_force"]
        np.testing.assert_allclose(force, source_force, rtol=1e-10, atol=1e-5)


def test_imported_heading_reaches_both_bodies_in_public_coupled_run(tmp_path):
    """The body-two source force survives the public shared-HDF5 path."""
    time = np.arange(201) * .1
    eta = .4 * np.cos(2 * np.pi * time / 7) + .15 * np.sin(2 * np.pi * time / 3)
    record = tmp_path / "etaData.mat"
    savemat(record, {"etaData": np.column_stack((time, eta))})

    wec = WEC("OSWEC two-body imported heading", body_to_body=True)
    flap = wec.body("flap", ROOT / "examples/data/oswec.h5",
                    mass=127_000, inertia=(1e7, 1e7, 1e7))
    wec.body("base", ROOT / "examples/data/oswec.h5")
    wec.coordinate("yaw", flap.move("yaw"))
    settings = dict(dt=.1, end_time=20, ramp_time=2, radiation_memory=10)
    oblique = wec.run(ImportedElevationWave(record, direction=15), **settings)
    head_on = wec.run(ImportedElevationWave(record, direction=0), **settings)

    for number in (1, 2):
        source_force = _imported_heading_source(number)[1]
        force = dict(oblique.raw.extra_outputs)[f"body{number}_excitation_force"]
        np.testing.assert_allclose(force, source_force, rtol=1e-10, atol=1e-5)
    np.testing.assert_array_equal(oblique.wave_elevation, head_on.wave_elevation)
    np.testing.assert_array_equal(oblique.bodies["base"].velocity, 0)
    # At oblique incidence the yaw excitation changes the independently
    # advanced response, while the imported elevation record stays fixed.
    assert np.max(np.abs(oblique.bodies["flap"].position[:, 5]
                         - head_on.bodies["flap"].position[:, 5])) > .01

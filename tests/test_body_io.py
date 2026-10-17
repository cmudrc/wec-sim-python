"""Read the RM3 hydrodynamic file included with the original port."""

from copy import deepcopy
from pathlib import Path
import sys

import numpy as np
import pytest
import scipy.io

ROOT = Path(__file__).resolve().parents[1]

from wecsim.hydrodynamics.bodyClass import BodyClass  # noqa: E402


def test_rm3_hydrodynamic_bodies_load():
    h5_file = ROOT / "examples/data/rm3.h5"
    for number, expected_name in ((1, "float"), (2, "spar")):
        body = BodyClass(str(h5_file))
        body.bodyNumber = number
        body.readH5file()
        assert body.name == expected_name
        assert body.hydroData["properties"]["name"] == expected_name
        assert body.hydroData["simulation_parameters"]["water_depth"] == "infinite"
        assert body.hydroData["hydro_coeffs"]["excitation"]["re"].shape[0] == 6


@pytest.mark.parametrize("source", ("hdf5", "memory"))
def test_rm3_regular_cic_force_preprocessing_matches_matlab_fixture(source):
    fixture = ROOT / "tests/test_objects/test_bodyclass/testData"
    case = fixture / "body_1_test"
    body = BodyClass(str(fixture / "hydroData/rm3.h5"))
    body.bodyNumber = 1
    body.readH5file()
    if source == "memory":
        # MCR can hand preloaded BEM data to the body instead of an HDF5 path.
        loaded = BodyClass(None)
        loaded.bodyNumber = 1
        loaded.loadHydroData(deepcopy(body.hydroData))
        body = loaded
    body.hydroStiffness = np.zeros((6, 6))
    body.viscDrag = {
        "Drag": np.zeros((6, 6)),
        "cd": np.zeros(6),
        "characteristicArea": np.zeros(6),
    }
    body.linearDamping = np.zeros((6, 6))
    body.mass = "equilibrium"
    body.hydroForcePre(
        0.785398163397448, [0], 601,
        np.loadtxt(case / "CTTime.txt").T, [], 0.1, 1000, 9.81,
        "regularCIC", np.loadtxt(case / "waveAmpTime.txt").T,
        1, [], 0, 0, 0,
    )

    for field, filename in (
        ("linearHydroRestCoef", "linearHydroRestCoef.txt"),
        ("fAddedMass", "fAddedMass.txt"),
    ):
        np.testing.assert_allclose(
            body.hydroForce[field], np.loadtxt(case / filename), rtol=1e-12,
        )
    np.testing.assert_allclose(
        body.hydroForce["fExt"]["re"],
        [18508.7581918053, 0.223532780361825, 1444829.62495620,
         6.41939080869329, 140467.443724649, 0.0412557412840770],
        rtol=1e-12,
    )
    np.testing.assert_allclose(
        body.hydroForce["fExt"]["im"],
        [546623.415011964, -0.0428753721889769, 447276.984830190,
         0.273118282507816, 4146814.04776144, 0.0138224113420678],
        rtol=1e-12,
    )
    np.testing.assert_allclose(
        body.hydroForce["irkb"], scipy.io.loadmat(case / "irkb.mat")["a"],
        rtol=1e-12, atol=1e-9,
    )


def test_oswec_directional_irregular_preprocessing_matches_matlab_fixture():
    fixture = ROOT / "tests/test_objects/test_bodyclass/testData"
    case = fixture / "body_3_test"

    def mat(name):
        values = scipy.io.loadmat(case / f"{name}.mat")
        return next(value for key, value in values.items() if not key.startswith("__"))

    body = BodyClass(str(fixture / "hydroData/oswec.h5"))
    body.bodyNumber = 1
    body.readH5file()
    body.hydroStiffness = np.zeros((6, 6))
    body.viscDrag = {
        "Drag": np.zeros((6, 6)),
        "cd": np.zeros(6),
        "characteristicArea": np.zeros(6),
    }
    body.linearDamping = np.zeros((6, 6))
    body.mass = 127000
    body.hydroForcePre(
        mat("w").T, [0, 30, 90], 301, mat("CTTime")[0], 500,
        0.1, 1000, 9.81, "irregular", mat("waveAmpTime").T,
        1, [], 0, 0, 0,
    )

    for field, expected in (
        ("linearHydroRestCoef", mat("linearHydroRestCoef")),
        ("fAddedMass", mat("fAddedMass")),
        ("irkb", mat("irkb")),
    ):
        np.testing.assert_allclose(body.hydroForce[field], expected, rtol=1e-10, atol=1e-8)
    for component in ("re", "im", "md"):
        np.testing.assert_allclose(
            body.hydroForce["fExt"][component], mat(component), rtol=1e-10, atol=1e-8,
        )


def test_rm3_irregular_preprocessing_matches_matlab_fixture():
    fixture = ROOT / "tests/test_objects/test_bodyclass/testData"
    case = fixture / "body_2_test"

    def mat(name):
        values = scipy.io.loadmat(case / f"{name}.mat")
        return next(value for key, value in values.items() if not key.startswith("__"))

    body = BodyClass(str(fixture / "hydroData/rm3.h5"))
    body.bodyNumber = 1
    body.readH5file()
    body.hydroStiffness = np.zeros((6, 6))
    body.viscDrag = {
        "Drag": np.zeros((6, 6)),
        "cd": np.zeros(6),
        "characteristicArea": np.zeros(6),
    }
    body.linearDamping = np.zeros((6, 6))
    body.mass = "equilibrium"
    body.hydroForcePre(
        np.loadtxt(case / "w.txt"), [0], 201,
        np.loadtxt(case / "CTTime.txt"), 500, 0.1, 1000, 9.81,
        "irregular", np.loadtxt(case / "waveAmpTime.txt").T,
        1, [], 0, 0, 0,
    )

    for field, filename in (
        ("linearHydroRestCoef", "linearHydroRestCoef.txt"),
        ("fAddedMass", "fAddedMass.txt"),
    ):
        np.testing.assert_allclose(
            body.hydroForce[field], np.loadtxt(case / filename), rtol=1e-10, atol=1e-8,
        )
    np.testing.assert_allclose(body.hydroForce["irkb"], mat("irkb"), rtol=1e-10, atol=1e-8)
    for component in ("re", "im", "md"):
        np.testing.assert_allclose(
            body.hydroForce["fExt"][component], mat(component), rtol=1e-10, atol=1e-8,
        )


@pytest.mark.parametrize("case_number,b2b", [(4, 0), (5, 1)])
@pytest.mark.parametrize("body_number", [1, 2])
def test_rm3_body_interaction_regular_preprocessing(case_number, b2b, body_number):
    fixture = ROOT / "tests/test_objects/test_bodyclass/testData"
    case = fixture / f"body_{case_number}_test"

    def mat(name):
        values = scipy.io.loadmat(case / f"{name}.mat")
        return next(value for key, value in values.items() if not key.startswith("__"))

    body = BodyClass(str(fixture / "hydroData/rm3.h5"))
    body.bodyNumber = body_number
    body.bodyTotal = [2]
    body.readH5file()
    body.hydroStiffness = np.zeros((6, 6))
    body.viscDrag = {
        "Drag": np.zeros((6, 6)),
        "cd": np.zeros(6),
        "characteristicArea": np.zeros(6),
    }
    body.linearDamping = np.zeros((6, 6))
    body.mass = "equilibrium"
    body.hydroForcePre(
        mat("w").T, [0], 601, mat("CTTime")[0], [],
        0.1, 1000, 9.81, "regular", mat("waveAmpTime").T,
        body_number, [], 0, 0, b2b,
    )
    prefix = f"body{body_number}_"
    for field in ("linearHydroRestCoef", "fAddedMass", "fDamping", "totDOF"):
        np.testing.assert_allclose(
            body.hydroForce[field], mat(prefix + field), rtol=1e-10, atol=1e-8,
        )
    for component in ("re", "im", "md"):
        np.testing.assert_allclose(
            body.hydroForce["fExt"][component], mat(prefix + component).ravel(),
            rtol=1e-10, atol=1e-8,
        )


@pytest.mark.parametrize("case_number,b2b,ss_calc", [
    (6, 0, 0), (7, 1, 0), (8, 0, 1), (9, 1, 1),
])
@pytest.mark.parametrize("body_number", [1, 2])
def test_rm3_body_interaction_cic_preprocessing(case_number, b2b, ss_calc, body_number):
    fixture = ROOT / "tests/test_objects/test_bodyclass/testData"
    case = fixture / f"body_{case_number}_test"

    def mat(name):
        values = scipy.io.loadmat(case / f"{name}.mat")
        return next(value for key, value in values.items() if not key.startswith("__"))

    body = BodyClass(str(fixture / "hydroData/rm3.h5"))
    body.bodyNumber = body_number
    body.bodyTotal = [2]
    body.readH5file()
    body.hydroStiffness = np.zeros((6, 6))
    body.viscDrag = {
        "Drag": np.zeros((6, 6)),
        "cd": np.zeros(6),
        "characteristicArea": np.zeros(6),
    }
    body.linearDamping = np.zeros((6, 6))
    body.mass = "equilibrium"
    body.hydroForcePre(
        mat("w").T, [0], 601, mat("CTTime")[0], [],
        0.1, 1000, 9.81, "regularCIC", mat("waveAmpTime").T,
        body_number, [], ss_calc, 0, b2b,
    )
    prefix = f"body{body_number}_"
    for field in ("linearHydroRestCoef", "fAddedMass", "irkb"):
        np.testing.assert_allclose(
            body.hydroForce[field], mat(prefix + field), rtol=1e-10, atol=1e-8,
        )
    for component in ("re", "im", "md"):
        np.testing.assert_allclose(
            body.hydroForce["fExt"][component], mat(prefix + component).ravel(),
            rtol=1e-10, atol=1e-8,
        )
    if ss_calc:
        for matrix in ("A", "B", "C", "D"):
            np.testing.assert_allclose(
                body.hydroForce["ssRadf"][matrix], mat(prefix + matrix),
                rtol=1e-10, atol=1e-8,
            )

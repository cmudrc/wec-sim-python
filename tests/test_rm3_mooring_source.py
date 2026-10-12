"""Check the pinned RM3 MooringMatrix source before pairing a Python solver."""

import os
from pathlib import Path

import numpy as np
import pytest
from scipy.io import loadmat

from wecsim.hydrodynamics.bodyClass import BodyClass


REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
APPLICATIONS = os.environ.get("WEC_SIM_APPLICATIONS_DIR")
pytestmark = pytest.mark.skipif(
    not (REFERENCE and APPLICATIONS),
    reason="paired MATLAB mooring output and Applications checkout not provided",
)


def _read(name):
    values = np.loadtxt(Path(REFERENCE) / name, delimiter=",")
    assert np.isfinite(values).all(), name
    return values


def test_mooring_matrix_source_exports_complete_motion_and_force():
    prefix = "RM3_MOORING_MATRIX"
    wave = _read(f"{prefix}_wave.csv")
    mooring = _read(f"{prefix}_mooring1.csv")
    stiffness = _read(f"{prefix}_stiffness.csv")
    assert stiffness.shape == (6, 6)
    np.testing.assert_allclose(stiffness, np.diag([100_000, 0, 0, 0, 0, 0]))
    assert mooring.shape == (40_001, 19)
    np.testing.assert_allclose(mooring[:, 0], np.arange(40_001) * 0.01,
                               rtol=0, atol=1e-8)
    assert wave.shape[1] == 2
    assert wave[0, 0] == 0
    assert wave[-1, 0] >= 400
    assert np.max(np.abs(wave[:, 1])) > 0.01
    assert np.max(np.abs(mooring[:, 13])) > 1_000
    for body in (1, 2):
        values = _read(f"{prefix}_MooringMatrix_body{body}.csv")
        assert values.shape == (40_001, 25)
        np.testing.assert_allclose(values[:, 0], mooring[:, 0],
                                   rtol=0, atol=1e-8)
    pto = _read(f"{prefix}_MooringMatrix_pto1.csv")
    assert pto.shape == (40_001, 25)
    np.testing.assert_allclose(pto[:, 0], mooring[:, 0], rtol=0, atol=1e-8)


def test_mooring_spring_and_imported_excitation_match_source():
    prefix = "RM3_MOORING_MATRIX"
    mooring = _read(f"{prefix}_mooring1.csv")
    stiffness = _read(f"{prefix}_stiffness.csv")
    expected_force = -mooring[:, 1:7] @ stiffness.T
    np.testing.assert_allclose(mooring[:, 13:19], expected_force,
                               rtol=0, atol=1e-5)

    apps = Path(APPLICATIONS)
    input_wave = loadmat(apps / "Mooring/MooringMatrix/etaData.mat")["etaData"]
    wave = _read(f"{prefix}_wave.csv")
    time = wave[:, 0]
    ramp = np.ones_like(time)
    early = time < 40
    ramp[early] = (1 - np.cos(np.pi * time[early] / 40)) / 2
    interpolated = np.interp(time, input_wave[:, 0], input_wave[:, 1])
    np.testing.assert_allclose(wave[:, 1], interpolated * ramp,
                               rtol=0, atol=1e-10)

    hydro = apps / "_Common_Input_Files/RM3/hydroData/rm3.h5"
    for number in (1, 2):
        body = BodyClass(str(hydro))
        body.bodyNumber = number
        body.bodyTotal = 2
        body.readH5file()
        body.hydroForce["userDefinedFe"] = np.zeros((len(time), 6))
        body.userDefinedExcitation(wave.T, 0.01, [0], 1000, 9.81)
        source = _read(f"{prefix}_MooringMatrix_body{number}.csv")
        # MATLAB ramps the imported elevation before convolution, then ramps
        # the resulting force once more when the body block applies it.
        calculated = body.hydroForce["userDefinedFe"] * ramp[:, None]
        assert np.max(np.abs(calculated - source[:, 19:25])) < 1e-3

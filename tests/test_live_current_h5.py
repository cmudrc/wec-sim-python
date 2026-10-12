"""Preprocess the HDF5 inputs in the pinned current MATLAB reference models."""

import os
from pathlib import Path
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]

from wecsim.hydrodynamics.bodyClass import BodyClass  # noqa: E402

MODEL = os.environ.get("WEC_SIM_REFERENCE_MODEL")
CORE = os.environ.get("WEC_SIM_MATLAB_CORE_DIR")
pytestmark = pytest.mark.skipif(
    not (MODEL in ("RM3", "OSWEC") and CORE),
    reason="current MATLAB reference-model HDF5 not provided",
)


@pytest.mark.parametrize("body_number", [1, 2])
def test_current_model_hydrodynamics_preprocess(body_number):
    h5_file = Path(CORE) / "examples" / MODEL / "hydroData" / f"{MODEL.lower()}.h5"
    body = BodyClass(str(h5_file))
    body.bodyNumber = body_number
    body.readH5file()
    expected_names = {"RM3": ("float", "spar"), "OSWEC": ("flap", "base")}
    assert body.name == expected_names[MODEL][body_number - 1]

    body.hydroStiffness = np.zeros((6, 6))
    body.viscDrag = {
        "Drag": np.zeros((6, 6)),
        "cd": np.zeros(6),
        "characteristicArea": np.zeros(6),
    }
    body.linearDamping = np.zeros((6, 6))
    body.mass = "equilibrium"
    if MODEL == "RM3":
        time = np.arange(601) * 0.1
        frequency = 2 * np.pi / 8
        directions = [0]
        num_frequencies = []
        wave_type = "regularCIC"
    else:
        time = np.arange(301) * 0.1
        frequency = np.linspace(0.0401, 19.999, 500)
        directions = [0, 30, 90]
        num_frequencies = 500
        wave_type = "irregular"
    body.hydroForcePre(
        frequency, directions, len(time), time, num_frequencies,
        0.1, 1000, 9.81, wave_type, np.vstack((time, np.zeros_like(time))),
        body_number, [], 0, 0, 0,
    )
    assert body.hydroForce["linearHydroRestCoef"].shape == (6, 6)
    assert body.hydroForce["fAddedMass"].shape == (6, 6)
    assert body.hydroForce["irkb"].shape == (len(time), 6, 6)
    for field in ("linearHydroRestCoef", "fAddedMass", "irkb"):
        assert np.isfinite(body.hydroForce[field]).all()
    assert np.isfinite(body.hydroForce["fExt"]["re"]).all()

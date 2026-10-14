"""Full-trace permeate flow through the published OSWEC RO membrane."""

import os
from pathlib import Path

import numpy as np
from wecsim import ReverseOsmosisMembrane


REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
FIXTURE = Path(__file__).parent / "fixtures/matlab/oswec_desalination/source.npz"


def test_published_membrane_permeate_flow():
    if REFERENCE:
        measured = np.loadtxt(
            Path(REFERENCE) / "OSWEC_DESALINATION_SOURCE_simout1.csv",
            delimiter=",",
        )
    else:
        with np.load(FIXTURE) as source:
            measured = source["simout1"]
    assert measured.shape == (30_001, 8)
    membrane = ReverseOsmosisMembrane(
        resistance=.6023e8,
        osmotic_pressure=30e5,
        regulation_range=5e4,
        max_area=.3,
        leakage_area=1e-12,
        discharge_coefficient=.7,
        fluid_density=850,
        laminar_pressure_ratio=.999,
    )
    actual = membrane.permeate_flow(measured[:, 6])
    np.testing.assert_allclose(actual, measured[:, 2], rtol=0, atol=1e-10)
    assert np.max(actual) > .04

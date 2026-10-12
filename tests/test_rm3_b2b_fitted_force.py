"""Replay source radiation force in the published RM3 Cases 5 and 6."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim.hydrodynamics.radiation import replay_rm3_fitted_radiation


APPLICATIONS = os.environ.get("WEC_SIM_APPLICATIONS_DIR")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
pytestmark = pytest.mark.skipif(
    not (APPLICATIONS and REFERENCE),
    reason="paired MATLAB RM3 state-space output not provided",
)


@pytest.mark.parametrize(("case", "body_to_body"), [(5, False), (6, True)])
def test_published_fitted_radiation_force(case, body_to_body):
    hydro = (Path(APPLICATIONS)
             / "_Common_Input_Files/RM3/hydroData/rm3.h5").resolve()
    traces = [np.loadtxt(
        Path(REFERENCE)
        / f"RM3_B2B_STATE_SPACE_B2B_Case{case}_body{body}.csv",
        delimiter=",",
    ) for body in (1, 2)]
    assert all(trace.shape == (4001, 37) for trace in traces)
    time = traces[0][:, 0]
    np.testing.assert_allclose(time, np.arange(4001) * 0.1,
                               rtol=0, atol=1e-10)
    np.testing.assert_allclose(traces[1][:, 0], time, rtol=0, atol=1e-10)
    velocity = np.stack([trace[:, 7:13] for trace in traces], axis=1)
    expected = np.stack([trace[:, 25:31] for trace in traces], axis=1)
    calculated = replay_rm3_fitted_radiation(
        hydro, time, velocity, body_to_body=body_to_body,
    )
    error = np.max(np.abs(calculated - expected), axis=0)
    peak = np.max(np.abs(expected), axis=0)
    assert np.all(error < 7e-4 * np.maximum(peak, 3)), (error, peak)

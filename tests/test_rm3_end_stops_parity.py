"""Check the effective source settings and forces of published RM3 End_Stops."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim.dynamics.hardStops import LinearHardStops


REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
pytestmark = pytest.mark.skipif(
    not REFERENCE, reason="paired MATLAB RM3 end-stop output not provided",
)


def test_pinned_end_stop_force_matches_logged_matlab_motion():
    reference = Path(REFERENCE)
    settings = np.loadtxt(reference / "RM3_END_STOPS_settings.csv", delimiter=",")
    np.testing.assert_allclose(
        settings, [-0.6, 0.6, 1e8, 1e8, 0, 0, 1e-4, 1e-4],
        rtol=0, atol=0,
    )
    pto = np.loadtxt(
        reference / "RM3_END_STOPS_End_Stops_pto1.csv", delimiter=",",
    )
    assert pto.shape == (4001, 49)
    stroke, speed, force = pto[:, 3], pto[:, 9], pto[:, 15]
    expected = (-1_200_000 * speed
                - 1e8 * np.maximum(stroke - 0.6, 0)
                + 1e8 * np.maximum(-0.6 - stroke, 0))
    assert np.max(np.abs(force - expected)) < 1e-5
    stops = LinearHardStops(-0.6, 0.6, 1e8, 1e8)
    assert np.max(np.abs(force - (-1_200_000 * speed
                                  + stops.force(stroke, speed)))) < 1e-5
    assert np.count_nonzero(stroke > 0.6) > 100
    assert np.count_nonzero(stroke < -0.6) > 100
    np.testing.assert_allclose(pto[:, 27], force, rtol=0, atol=1e-8)

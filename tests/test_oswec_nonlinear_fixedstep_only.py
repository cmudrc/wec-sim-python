"""Measure the source motion when only Simulink ode4's step is refined."""

import os
from pathlib import Path

import numpy as np
import pytest
from scipy.io import loadmat


FIXTURE = (Path(__file__).parent / "fixtures" / "matlab"
           / "oswec_nonlinear_viz_published_motion.npz")


@pytest.mark.skipif(
    not os.getenv("WEC_SIM_OSWEC_FIXEDSTEP_ONLY_DIR"),
    reason="pinned OSWEC fixed-step-only source output not provided",
)
def test_published_oswec_with_only_fixed_solver_step_refined():
    result = loadmat(Path(os.environ["WEC_SIM_OSWEC_FIXEDSTEP_ONLY_DIR"])
                     / "solver-only.mat")
    with np.load(FIXTURE) as source:
        published = {key: source[key] for key in source.files}

    for key in ("dt", "dtOut", "cicDt", "nonlinearDt"):
        assert float(result[key].item()) == .1
    assert float(result["fixedStep"].item()) == .0125
    np.testing.assert_allclose(result["time"].ravel(),
                               published["time"].ravel(),
                               rtol=0, atol=1e-10)
    for key in ("pose", "velocity", "acceleration", "forceExcitation",
                "forceRadiationDamping", "forceAddedMass", "forceRestoring"):
        assert result[key].shape == published[key].shape == (1201, 6)
        assert np.isfinite(result[key]).all()

    # The new arm is a controlled source diagnostic, not a Python parity
    # gate. Report its response to the solver-only change without assuming
    # which internal source setting dominates beforehand.
    for label, key, column in (
        ("pitch", "pose", 4), ("surge", "pose", 0),
        ("heave", "pose", 2), ("pitch speed", "velocity", 4),
        ("restoring pitch moment", "forceRestoring", 4),
        ("excitation pitch moment", "forceExcitation", 4),
        ("added-mass pitch moment", "forceAddedMass", 4),
    ):
        difference = result[key][:, column] - published[key][:, column]
        print(f"{label}: published-to-solver-only max="
              f"{np.max(np.abs(difference)):.9g}, "
              f"rms={np.sqrt(np.mean(difference**2)):.9g}")

"""Portable independent ellipsoid ode45 gate against checked R2025b source."""

from pathlib import Path

import numpy as np
import pytest

from test_ellipsoid_ode45_convergence import (
    OUTPUT_STEPS, STEPS, _validate_full_gate,
)


FIXTURE = (Path(__file__).resolve().parents[1] / "tests/fixtures/matlab"
           / "ellipsoid_ode45_convergence")


@pytest.mark.parametrize("case", ["Regular", "RegularCIC"])
def test_archived_source_screen_and_full_duration_independent_motion(case):
    assert (FIXTURE / "source.npz").is_file()
    assert (FIXTURE / "hydroData/ellipsoid.h5").is_file()
    assert (FIXTURE / "geometry/elipsoid.stl").is_file()
    key = "regular_cic" if case == "RegularCIC" else "regular"
    with np.load(FIXTURE / "source.npz", allow_pickle=False) as archived:
        for label, step in STEPS.items():
            np.testing.assert_allclose(
                archived[f"{key}_{label}_settings"],
                [[step, OUTPUT_STEPS[label], step, step, step]],
                rtol=0, atol=1e-12,
            )
        source = {label: archived[f"{key}_{label}"] for label in STEPS}
    _validate_full_gate(source, case, FIXTURE)

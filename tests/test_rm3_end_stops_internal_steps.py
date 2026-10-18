"""Locate the published End_Stops delay residual on logged solver steps."""

import os
from pathlib import Path

import numpy as np
import pytest


FIXTURE = Path(__file__).parent / "fixtures/matlab/rm3_end_stops_internal_steps.npz"
BASELINE = Path(__file__).parent / "fixtures/matlab/rm3_end_stops_maxstep_400s.npz"
OUTPUT = os.environ.get("WEC_SIM_RM3_END_STOPS_MAXSTEP_DIR")


def _check(data, *, pinned):
    grid = data["output_time"]
    internal = data["internal_times"]
    np.testing.assert_allclose(grid, np.arange(4001) * .1, rtol=0, atol=1e-9)
    assert np.isfinite(internal).all() and np.all(np.diff(internal) > 0)
    assert internal[0] == 0 and internal[-1] == 400
    # Every published output sample is a logged solver time. Some solver
    # steps after stop contact are absent from the output histories.
    indices = np.searchsorted(internal, grid - 1e-9)
    np.testing.assert_allclose(internal[indices], grid, rtol=0, atol=1e-9)
    off_grid = internal[np.abs(internal - np.round(internal * 10) / 10) > 1e-8]
    assert len(off_grid) > 0

    # A saved force at t uses acceleration history around t - 100 ns.
    # Two saved accelerations suffice only if no extra solver step occurred
    # within the preceding two output intervals.
    lower = np.searchsorted(off_grid, grid[:-2], side="right")
    upper = np.searchsorted(off_grid, grid[2:], side="left")
    impacted = upper > lower
    for body in (1, 2):
        acceleration = data[f"body{body}_acceleration"]
        applied_force = data[f"body{body}_added_force"]
        applied_mass = data[f"body{body}_added_matrix"]
        assert acceleration.shape == applied_force.shape == (4001, 6)
        assert applied_mass.shape == (6, 6)
        assert all(np.isfinite(value).all() for value in
                   (acceleration, applied_force, applied_mass))
        delayed = acceleration[1:-1] + (1 - 1e-7 / .1) * (
            acceleration[1:-1] - acceleration[:-2]
        )
        predicted = delayed @ applied_mass.T
        residual = np.max(np.abs(predicted[:, [0, 2]]
                                 - applied_force[2:, [0, 2]]), axis=1)
        assert residual[~impacted].max() < .001
        assert not np.any((residual > .001) & ~impacted)
        if pinned:
            assert len(off_grid) == 81
            assert np.count_nonzero(impacted) == 78
            np.testing.assert_array_equal(residual > .001, impacted)
            assert residual.max() > 17_000_000


def test_archived_internal_step_history_explains_force_residual():
    with np.load(FIXTURE, allow_pickle=False) as data, \
         np.load(BASELINE, allow_pickle=False) as baseline:
        np.testing.assert_array_equal(data["published_pto"],
                                      baseline["published_pto"])
        _check(data, pinned=True)


@pytest.mark.skipif(not OUTPUT, reason="fresh MATLAB MaxStep audit absent")
def test_fresh_internal_step_history_explains_force_residual():
    directory = Path(OUTPUT)
    pto = np.loadtxt(directory / "published_pto.csv", delimiter=",")
    assert pto.shape == (4001, 4)
    with np.load(BASELINE, allow_pickle=False) as baseline:
        for column, limit in enumerate((1e-7, 1e-6, 1e-2)):
            np.testing.assert_allclose(
                pto[:, column + 1], baseline["published_pto"][:, column],
                rtol=0, atol=limit,
            )
    arrays = {
        "output_time": pto[:, 0],
        "internal_times": np.loadtxt(directory / "published_internal_times.csv",
                                     delimiter=","),
    }
    for body in (1, 2):
        for key, suffix in (("acceleration", "acceleration"),
                            ("added_force", "added_force"),
                            ("added_matrix", "added_matrix")):
            arrays[f"body{body}_{key}"] = np.loadtxt(
                directory / f"published_body{body}_{suffix}.csv", delimiter=",",
            )
    assert arrays["output_time"].shape == (4001,)
    _check(arrays, pinned=False)

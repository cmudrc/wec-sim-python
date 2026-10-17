"""Verify that a logged MoorDyn callback preserves the source trajectory."""

import os
from pathlib import Path

import numpy as np
import pytest


OUTPUT = os.getenv("WEC_SIM_RM3_CALLBACK_DIR")
FIXTURE = (Path(__file__).parent / "fixtures/matlab"
           / "rm3_moordyn_initial_step_2s/signals.npz")


@pytest.mark.skipif(not OUTPUT, reason="MATLAB callback audit not provided")
def test_logged_callback_does_not_change_source_motion():
    directory = Path(OUTPUT)
    callback_traces = {}
    with np.load(FIXTURE, allow_pickle=False) as expected:
        for label, prefix in (("auto", "default"),
                              ("initial_001", "initial_step")):
            phase = np.loadtxt(directory / f"{label}_phase.csv", delimiter=",")
            wave = np.loadtxt(directory / f"{label}_wave.csv", delimiter=",")
            np.testing.assert_array_equal(phase, expected[f"{prefix}_phase"])
            np.testing.assert_allclose(wave[:, 0], expected["time"],
                                       rtol=0, atol=1e-12)
            np.testing.assert_allclose(wave[:, 1], expected[f"{prefix}_wave"],
                                       rtol=0, atol=1e-12)
            for body in (0, 1):
                trace = np.loadtxt(directory / f"{label}_body{body + 1}.csv",
                                   delimiter=",")
                assert trace.shape == (201, 13)
                for name, column, limit in (("position", 1, 1e-4),
                                            ("velocity", 7, 1e-3)):
                    difference = np.max(np.abs(
                        trace[:, column:column + 6]
                        - expected[f"{prefix}_{name}"][body]
                    ))
                    assert difference < limit, (label, body, name, difference)

            callbacks = np.loadtxt(
                directory / f"{label}_moordyn_calls.csv", delimiter=",")
            assert callbacks.ndim == 2 and callbacks.shape[1] == 20
            assert np.isfinite(callbacks).all()
            assert np.all(callbacks[:, 1] >= 0)
            callback_traces[label] = callbacks
            print(f"{label}: {len(callbacks)} MoorDyn calls; "
                  f"first event={sum(callbacks[0, :2]):.9g} s, "
                  f"last event={sum(callbacks[-1, :2]):.9g} s")

    default = callback_traces["auto"]
    fixed = callback_traces["initial_001"]
    assert default.shape == fixed.shape == (200, 20)
    # The second MATLAB arm changes ode45's startup step, not MoorDyn's
    # coupling interval or the number/times of its native calls.
    for calls in (default, fixed):
        np.testing.assert_allclose(calls[:, 1], .01, rtol=0, atol=1e-12)
        np.testing.assert_allclose(calls[:, 0] + calls[:, 1],
                                   np.arange(1, 201) * .01,
                                   rtol=0, atol=1e-12)
    np.testing.assert_allclose(default[:, :2], fixed[:, :2],
                               rtol=0, atol=1e-12)

    # At the first call, connection state already differs. This is an
    # upstream solver-state difference; the returned load then amplifies it.
    first_pose = abs(default[0, 2] - fixed[0, 2])
    first_speed = abs(default[0, 8] - fixed[0, 8])
    first_moment = abs(default[0, 18] - fixed[0, 18])
    assert 1e-5 < first_pose < 2e-5
    assert .003 < first_speed < .004
    assert 10_000 < first_moment < 20_000

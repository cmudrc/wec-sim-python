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
            print(f"{label}: {len(callbacks)} MoorDyn calls; "
                  f"first time={callbacks[0, 0]:.9g} s, "
                  f"last time={callbacks[-1, 0]:.9g} s")

"""The full WaveStar controller does not retain a 50 ms command grid at 25 ms ode8."""

from pathlib import Path

import numpy as np


FIXTURE = Path(__file__).parent / "fixtures/matlab/wavestar_nmpc_full_step"
PRECONTROL = Path(__file__).parent / "fixtures/matlab/wavestar_fixed_input_step"


def test_archived_published_input_matches_the_precontrol_fixture():
    with np.load(FIXTURE / "baseline.npz") as published, \
            np.load(PRECONTROL / "source.npz") as short:
        assert published["body"].shape == (4501, 25)
        assert published["pto"].shape == (4501, 43)
        assert published["components"].shape == (500, 4)
        np.testing.assert_array_equal(
            published["body"][:200],
            short["50ms_WECCCOMP_Nonlinear_Model_Predictive_body1"],
        )
        np.testing.assert_array_equal(
            published["pto"][:200],
            short["50ms_WECCCOMP_Nonlinear_Model_Predictive_pto1"],
        )
        np.testing.assert_array_equal(
            published["components"], short["50ms_components"],
        )
        np.testing.assert_array_equal(published["wave"][:200], short["50ms_wave"])


def test_full_25ms_solver_run_changes_the_logged_controller_cadence():
    with np.load(FIXTURE / "baseline.npz") as published, \
            np.load(FIXTURE / "solver_25ms.npz") as refined:
        # These inputs are identical; the source time grid is the variable.
        np.testing.assert_array_equal(refined["wave"], published["wave"])
        np.testing.assert_array_equal(refined["components"],
                                      published["components"])
        assert refined["command"].shape == (9001,)
        for signal in ("command", "power", "energy"):
            np.testing.assert_allclose(
                published[f"{signal}_time"], np.arange(4501) * .05,
                rtol=0, atol=1e-9,
            )
            np.testing.assert_allclose(
                refined[f"{signal}_time"], np.arange(9001) * .025,
                rtol=0, atol=1e-9,
            )

        # The command itself changes between old 50 ms ticks, so taking
        # every second logged sample cannot turn this into a solver-only audit.
        command = refined["command"]
        half_step_change = np.abs(command[1::2] - command[:-1:2])
        after_30s = refined["command_time"][1::2] >= 30
        assert np.count_nonzero(half_step_change[after_30s] > 1e-6) > 0.99 * 3900
        assert np.max(half_step_change[after_30s]) > 1.5
        assert np.max(np.abs(command[::2] - published["command"])) > 14
        assert abs(refined["energy"][-1] - published["energy"][-1]) > 30

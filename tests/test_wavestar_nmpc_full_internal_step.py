"""The full WaveStar controller does not retain a 50 ms command grid at 25 ms ode8."""

from pathlib import Path

import numpy as np

from wecsim import WaveStarNmpcPTO
from wecsim.dynamics.wavestar import run_wavestar_published
from wecsim.environment.irregularWave import jonswap_equal_energy_components


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


def test_full_source_motion_changes_with_controller_cadence():
    """The 25 ms source is a different coupled trajectory, not a plant control."""
    with np.load(FIXTURE / "baseline.npz") as published, \
            np.load(FIXTURE / "solver_25ms.npz") as refined, \
            np.load(PRECONTROL / "source.npz") as short:
        body, pto = refined["body"], refined["pto"]
        assert body.shape == (4501, 25)
        assert pto.shape == (9001, 43)
        np.testing.assert_allclose(body[:, 0], np.arange(4501) * .05,
                                   rtol=0, atol=1e-9)
        np.testing.assert_allclose(pto[:, 0], np.arange(9001) * .025,
                                   rtol=0, atol=1e-9)
        np.testing.assert_allclose(
            body[:200], short["25ms_WECCCOMP_Nonlinear_Model_Predictive_body1"],
            rtol=0, atol=1e-9,
        )
        np.testing.assert_allclose(
            pto[:399], short["25ms_WECCCOMP_Nonlinear_Model_Predictive_pto1"],
            rtol=0, atol=1e-9,
        )
        assert np.max(np.abs(body[:, 5] - published["body"][:, 5])) > .16
        assert np.max(np.abs(body[:, 11] - published["body"][:, 11])) > 1.18
        assert np.max(np.abs(pto[::2, 33] - published["pto"][:, 33])) > 70
        old_work = np.trapezoid(
            published["pto"][:, 33] * published["pto"][:, 9],
            x=published["pto"][:, 0],
        )
        new_work = np.trapezoid(pto[:, 33] * pto[:, 9], x=pto[:, 0])
        assert .20 < abs(new_work - old_work) / abs(old_work) < .21


def test_independent_python_is_closer_to_published_than_25ms_source():
    """A changed controller clock cannot be used as a refined parity target."""
    with np.load(FIXTURE / "baseline.npz") as published, \
            np.load(FIXTURE / "solver_25ms.npz") as refined:
        components = published["components"]
        np.testing.assert_array_equal(components, refined["components"])
        sea = jonswap_equal_energy_components(
            PRECONTROL / "wavestar.h5", significant_height=.1042,
            peak_period=1.836, directions=np.array([0.]),
            spreading=np.array([1.]), gamma=3.3,
            phase=components[:, 3, None],
        )
        controller = WaveStarNmpcPTO(plant_dt=.001, control_dt=.05)
        python = run_wavestar_published(
            PRECONTROL / "wavestar.h5", sea, dt=.001, end_time=225,
            ramp_time=25, g=9.80665, pto_controller=controller,
            output_stride=50,
        )
        np.testing.assert_allclose(python.time, published["body"][:, 0],
                                   rtol=0, atol=1e-9)

        def maximum(a, b):
            return np.max(np.abs(a - b))

        old_pitch_gap = maximum(python.angle, published["body"][:, 5])
        new_pitch_gap = maximum(python.angle, refined["body"][:, 5])
        old_force_gap = maximum(python.pto_force, published["pto"][:, 33])
        new_force_gap = maximum(python.pto_force, refined["pto"][::2, 33])
        assert old_pitch_gap < .05 < .15 < new_pitch_gap
        assert old_force_gap < 15 < 60 < new_force_gap
        python_work = np.trapezoid(python.pto_force * python.pto_speed,
                                   x=python.time)
        old_work = np.trapezoid(
            published["pto"][:, 33] * published["pto"][:, 9],
            x=published["pto"][:, 0],
        )
        pto = refined["pto"]
        new_work = np.trapezoid(pto[:, 33] * pto[:, 9], x=pto[:, 0])
        assert abs(python_work - old_work) / abs(old_work) < .07
        assert abs(python_work - new_work) / abs(new_work) > .15

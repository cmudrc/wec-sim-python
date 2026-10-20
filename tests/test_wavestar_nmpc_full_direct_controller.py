"""Gate a full-duration NMPC QP replay on saved MATLAB controller inputs."""

from pathlib import Path

import numpy as np

from wecsim import WaveStarNmpcController


SOURCE = (Path(__file__).parent / "fixtures/matlab/wavestar_nmpc_full_step"
          / "direct_controller_inputs.npz")
INDEPENDENT_INPUT = (Path(__file__).parent
                     / "fixtures/matlab/wavestar_nmpc_independent_qp/source.npz")
EARLY_HORIZONS = (Path(__file__).parent
                  / "fixtures/matlab/wavestar_nmpc_independent_qp/early_horizons.npz")


def test_full_duration_direct_nmpc_controller_replay():
    with np.load(SOURCE, allow_pickle=False) as saved:
        time = saved["time"]
        state = saved["estimated_state"]
        forecast = saved["forecast"]
        source = saved["source_command"]
        direct = saved["direct_command"]
    assert state.shape == (4201, 5) and forecast.shape == (40, 4201)
    np.testing.assert_allclose(time, 15 + np.arange(4201) * .05,
                               rtol=0, atol=1e-9)
    assert np.isfinite(state).all() and np.isfinite(forecast).all()
    assert np.max(np.abs(source)) > 11.9  # The replay reaches torque saturation.
    assert np.max(np.abs(direct - source)) < 1e-9

    controller = WaveStarNmpcController()
    python_command = np.array([
        controller.step(state[index], forecast[:, index], float(instant))
        for index, instant in enumerate(time)
    ])
    gap = np.abs(python_command - direct)
    assert np.max(gap[time <= 20]) < 5e-6
    assert np.max(gap) < 1e-3
    assert np.sqrt(np.mean(gap**2)) < 5e-5


def test_independent_motion_inputs_against_matlab_quadprog():
    """Pair QP solvers on states/forecasts from Python's own 30 s motion."""
    with np.load(INDEPENDENT_INPUT, allow_pickle=False) as saved:
        time = saved["time"]
        state = saved["estimated_state"]
        forecast = saved["forecast"]
        direct = saved["matlab_command"]
        original_python = saved["python_command"]
    assert state.shape == (4201, 5) and forecast.shape == (40, 4201)
    np.testing.assert_allclose(time, 15 + np.arange(4201) * .05,
                               rtol=0, atol=1e-12)
    assert np.isfinite(state).all() and np.isfinite(forecast).all()
    assert np.count_nonzero(np.abs(direct) > 11.9) >= 40  # Saturated QPs.

    controller = WaveStarNmpcController()
    python_command = np.array([
        controller.step(state[index], forecast[:, index], float(instant))
        for index, instant in enumerate(time)
    ])
    # The replay was generated from the same Python states and forecasts;
    # MATLAB's pinned NMPC.m recomputed every command with quadprog.
    np.testing.assert_allclose(python_command, original_python,
                               rtol=0, atol=1e-3)
    gap = np.abs(python_command - direct)
    assert np.max(gap[time <= 20]) < 5e-5
    assert np.max(gap) < 1e-3
    assert np.sqrt(np.mean(gap**2)) < 5e-5


def test_early_independent_inputs_pair_complete_nmpc_plan():
    """Check all 40 planned torques where autonomous feedback first separates."""
    with np.load(INDEPENDENT_INPUT, allow_pickle=False) as saved:
        time = saved["time"][:40]
        state = saved["estimated_state"][:40]
        forecast = saved["forecast"][:, :40]
        matlab_command = saved["matlab_command"][:40]
    with np.load(EARLY_HORIZONS, allow_pickle=False) as saved:
        source_time = saved["time"]
        matlab_plan = saved["matlab_plan"]
    assert matlab_plan.shape == (40, 40)
    assert np.isfinite(matlab_plan).all()
    np.testing.assert_allclose(source_time, time, rtol=0, atol=1e-12)
    np.testing.assert_allclose(matlab_plan[:, 0], matlab_command,
                               rtol=0, atol=1e-9)

    controller = WaveStarNmpcController()
    python_plan = np.empty_like(matlab_plan)
    for index, instant in enumerate(time):
        controller.step(state[index], forecast[:, index], float(instant))
        python_plan[index] = controller._nominal
    # The existing early-window direct-QP gate uses 5e-5 N m. Apply that
    # same limit to every future command in each independently solved plan.
    assert np.max(np.abs(python_plan - matlab_plan)) < 5e-5

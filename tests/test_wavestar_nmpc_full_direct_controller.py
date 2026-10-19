"""Gate a full-duration NMPC QP replay on saved MATLAB controller inputs."""

from pathlib import Path

import numpy as np

from wecsim import WaveStarNmpcController


SOURCE = (Path(__file__).parent / "fixtures/matlab/wavestar_nmpc_full_step"
          / "direct_controller_inputs.npz")
INDEPENDENT_INPUT = (Path(__file__).parent
                     / "fixtures/matlab/wavestar_nmpc_independent_qp/source.npz")


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
    # MATLAB's pinned NMPC.m recomputed every command with quadprog on these
    # independently generated Python states and forecasts.
    gap = np.abs(python_command - direct)
    assert np.max(gap[time <= 20]) < 1e-5
    assert np.max(gap) < 2.5e-4
    assert np.sqrt(np.mean(gap**2)) < 1e-5

"""Gate a full-duration NMPC QP replay on saved MATLAB controller inputs."""

from pathlib import Path

import numpy as np

from wecsim import WaveStarNmpcController


SOURCE = (Path(__file__).parent / "fixtures/matlab/wavestar_nmpc_full_step"
          / "direct_controller_inputs.npz")


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

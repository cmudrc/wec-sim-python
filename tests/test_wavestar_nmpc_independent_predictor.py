"""Pair MATLAB AR forecasts with the same independently generated Python inputs."""

from pathlib import Path

import numpy as np
from scipy.io import loadmat

from wecsim import WaveStarNmpcObserver, WaveStarNmpcPredictor
from wecsim.dynamics.wavestar import WaveStarLinkage


FIXTURE = (Path(__file__).parent
           / "fixtures/matlab/wavestar_nmpc_independent_predictor")
QP_INPUT = (Path(__file__).parent
            / "fixtures/matlab/wavestar_nmpc_independent_qp/source.npz")


def test_full_duration_matlab_forecast_on_independent_python_motion():
    with np.load(FIXTURE / "python_motion.npz", allow_pickle=False) as motion:
        time = motion["time"]
        pitch = motion["pitch"]
        command = motion["command"]
    assert time.shape == pitch.shape == command.shape == (4501,)
    np.testing.assert_allclose(time, np.arange(4501) * .05,
                               rtol=0, atol=1e-10)
    assert np.isfinite(pitch).all() and np.isfinite(command).all()

    # Reconstruct the observer input from Python's saved autonomous motion.
    # The active states must also identify the independent QP replay fixture.
    observer = WaveStarNmpcObserver()
    stroke = WaveStarLinkage().pto_stroke(pitch)
    states = []
    previous = 0.
    for position, request in zip(stroke, command):
        states.append(observer.step(float(position), float(previous)))
        previous = request
    states = np.asarray(states)
    with np.load(QP_INPUT, allow_pickle=False) as qp:
        np.testing.assert_allclose(time[300:], qp["time"],
                                   rtol=0, atol=1e-10)
        np.testing.assert_allclose(command[300:], qp["python_command"],
                                   rtol=0, atol=1e-12)
        np.testing.assert_allclose(states[300:], qp["estimated_state"],
                                   rtol=0, atol=1e-12)

    matlab_input = np.loadtxt(FIXTURE / "python_moment.csv", delimiter=",")
    assert matlab_input.shape == (4501, 2)
    np.testing.assert_allclose(matlab_input[:, 0], time,
                               rtol=0, atol=1e-10)
    np.testing.assert_allclose(matlab_input[:, 1], states[:, -1],
                               rtol=0, atol=1e-12)

    predictor = WaveStarNmpcPredictor()
    forecast = np.array([predictor.step(float(moment))
                         for moment in states[:, -1]])
    refits = np.arange(200, 4501, 20)  # 10–225 s, every source refit.
    matlab = np.loadtxt(FIXTURE / "matlab_forecast.csv", delimiter=",")
    assert len(refits) == 216 and matlab.shape == (216, 41)
    np.testing.assert_allclose(matlab[:, 0], time[refits],
                               rtol=0, atol=1e-10)
    assert np.isfinite(matlab).all()
    assert np.max(np.abs(matlab[:, 1:])) > 10
    gap = forecast[refits] - matlab[:, 1:]
    assert np.max(np.abs(gap)) < 1e-6  # N m, all 8,640 forecasts.
    assert np.sqrt(np.mean(gap**2)) < 1e-7  # N m RMS.

    # MATLAB's pinned predictor.m holds its fitted AR model between one-second
    # refits. The controller consumes every 0.05 s forecast from 15 s onward.
    active = loadmat(FIXTURE / "matlab_active_forecast.mat")["activeRows"]
    assert active.shape == (4201, 41)
    np.testing.assert_allclose(active[:, 0], time[300:],
                               rtol=0, atol=1e-10)
    assert np.isfinite(active).all()
    active_gap = forecast[300:] - active[:, 1:]
    assert np.max(np.abs(active_gap)) < 1e-6  # N m, 168,040 values.
    assert np.sqrt(np.mean(active_gap**2)) < 1e-7  # N m RMS.

"""Pair the published WaveStar NMPC signal path without a MATLAB runner.

This is a component gate on saved stroke and preceding torque requests. It
does not establish independently coupled published-case motion parity.
"""

from pathlib import Path

import numpy as np

from wecsim import (
    WaveStarNmpcActuator,
    WaveStarNmpcController,
    WaveStarNmpcObserver,
    WaveStarNmpcPredictor,
)


FIXTURE = (Path(__file__).resolve().parent / "fixtures" / "matlab"
           / "wavestar_nmpc" / "controller_20s.npz")


def test_published_nmpc_observer_forecast_command_and_actuator_prefix():
    with np.load(FIXTURE) as source:
        time = source["time"]
        stroke = source["stroke"]
        command = source["command"]
        state = source["estimated_state"]
        forecast = source["forecast"]
        force = source["axial_force"]

    assert time.shape == stroke.shape == command.shape == force.shape == (401,)
    assert state.shape == (401, 5)
    assert forecast.shape == (40, 401)
    assert all(np.isfinite(values).all()
               for values in (time, stroke, command, state, forecast, force))
    np.testing.assert_allclose(time, np.arange(401) * .05,
                               rtol=0, atol=1e-12)
    assert np.max(np.abs(command[:200])) < 1e-12
    assert np.max(np.abs(command[300:])) > 1
    assert np.max(np.abs(forecast[:, 200:220])) > 1e6

    observer = WaveStarNmpcObserver()
    predictor = WaveStarNmpcPredictor()
    controller = WaveStarNmpcController()
    actuator = WaveStarNmpcActuator()
    actual_state = []
    actual_forecast = []
    actual_command = []
    actual_force = []
    for index, instant in enumerate(time):
        # The published preceding request is an input to the observer;
        # everything after that is independently recomputed in Python.
        estimated = observer.step(stroke[index], command[index - 1]
                                  if index else 0.)
        predicted = predictor.step(estimated[-1])
        actual_state.append(estimated)
        actual_forecast.append(predicted)
        actual_command.append(controller.step(estimated, predicted, instant))
        actual_force.append(actuator.step(command[index], stroke[index]))

    actual_state = np.asarray(actual_state)
    actual_forecast = np.asarray(actual_forecast).T
    actual_command = np.asarray(actual_command)
    actual_force = np.asarray(actual_force)
    assert np.max(np.abs(actual_state - state)) < 1e-10
    np.testing.assert_allclose(actual_forecast[:, :200], forecast[:, :200],
                               rtol=0, atol=1e-12)
    # The first AR fit is unstable at source scale; the matched source test
    # uses these same absolute limits for composed observer/predictor replay.
    assert np.max(np.abs(actual_forecast[:, 200:220]
                         - forecast[:, 200:220])) < 25
    assert np.max(np.abs(actual_forecast[:, 220:]
                         - forecast[:, 220:])) < 1e-9
    command_error = actual_command - command
    assert np.max(np.abs(command_error)) < .01
    assert np.sqrt(np.mean(command_error**2)) < .0003
    assert np.max(np.abs(actual_force - force)) < 1e-7

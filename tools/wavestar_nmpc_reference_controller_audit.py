"""Replay a checked WaveStar source controller log on its own stroke.

This source-driven component audit never supplies source motion or requests to
an autonomous Python plant. It isolates observer, AR, and QP replay errors.
"""

import argparse
import json
import platform
from hashlib import sha256
from pathlib import Path

import numpy as np
import scipy

from wecsim import (
    WaveStarNmpcController, WaveStarNmpcObserver, WaveStarNmpcPredictor,
)


COMMAND_MAX_LIMIT = .01
COMMAND_RMS_LIMIT = .0003


def _stats(values, time):
    difference = np.abs(np.asarray(values))
    return {
        "max_abs": float(difference.max()),
        "time_of_max_s": float(time[difference.argmax()]),
        "rms_sample_max": float(np.sqrt(np.mean(difference**2))),
        "max_active_after_15s": float(difference[300:].max()),
        "max_30_to_225s": float(difference[600:].max()),
    }


def audit(source_path, controller_path, expected_source_sha256=None,
          expected_controller_sha256=None):
    source_digest = sha256(source_path.read_bytes()).hexdigest()
    controller_digest = sha256(controller_path.read_bytes()).hexdigest()
    if (expected_source_sha256 is not None and
            source_digest != expected_source_sha256):
        raise ValueError("The source trace hash differs from the checked input")
    if (expected_controller_sha256 is not None and
            controller_digest != expected_controller_sha256):
        raise ValueError("The controller log hash differs from the checked input")
    with np.load(source_path, allow_pickle=False) as source:
        time = source["time"].copy()
        stroke = source["stroke"].copy()
        source_request = source["command"].copy()
    with np.load(controller_path, allow_pickle=False) as controller_log:
        log_time = controller_log["time"].copy()
        logged_state = controller_log["estimated_state"].copy()
        logged_forecast = controller_log["forecast"].copy()
        logged_request = controller_log["command"].copy()
    if (time.shape != (4501,) or stroke.shape != time.shape or
            source_request.shape != time.shape or log_time.shape != time.shape or
            logged_state.shape != (4501, 5) or
            logged_forecast.shape != (4501, 40) or
            logged_request.shape != time.shape):
        raise ValueError("The checked WaveStar controller log has wrong shapes")
    if not all(np.isfinite(array).all() for array in (
        time, stroke, source_request, log_time, logged_state,
        logged_forecast, logged_request,
    )):
        raise ValueError("The checked WaveStar controller log is nonfinite")
    np.testing.assert_allclose(time, np.arange(4501) * .05,
                               rtol=0, atol=1e-9)
    np.testing.assert_array_equal(log_time, time)
    np.testing.assert_array_equal(logged_request, source_request)

    observer = WaveStarNmpcObserver(dt=.05)
    predictor = WaveStarNmpcPredictor(dt=.05)
    complete_controller = WaveStarNmpcController(dt=.05)
    frozen_input_qp = WaveStarNmpcController(dt=.05)
    state_gap = np.empty(len(time))
    forecast_gap = np.empty(len(time))
    complete_gap = np.empty(len(time))
    qp_gap = np.empty(len(time))
    for index, instant in enumerate(time):
        prior = source_request[index - 1] if index else 0.
        state = observer.step(float(stroke[index]), float(prior))
        forecast = predictor.step(float(state[-1]))
        complete_request = complete_controller.step(state, forecast,
                                                    float(instant))
        qp_request = frozen_input_qp.step(logged_state[index],
                                          logged_forecast[index],
                                          float(instant))
        state_gap[index] = np.max(np.abs(state - logged_state[index]))
        forecast_gap[index] = np.max(np.abs(
            forecast - logged_forecast[index]
        ))
        complete_gap[index] = abs(complete_request - source_request[index])
        qp_gap[index] = abs(qp_request - source_request[index])

    complete = _stats(complete_gap, time)
    frozen = _stats(qp_gap, time)
    return {
        "source_trace_sha256": source_digest,
        "controller_signals_sha256": controller_digest,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "samples": len(time),
        "observer_state": _stats(state_gap, time),
        "ar_forecast": _stats(forecast_gap, time),
        "complete_controller_request": complete,
        "frozen_source_state_forecast_qp_request": frozen,
        "command_component_limits": {
            "max_abs_Nm": COMMAND_MAX_LIMIT,
            "rms_Nm": COMMAND_RMS_LIMIT,
        },
        "complete_controller_component_pass": (
            complete["max_abs"] < COMMAND_MAX_LIMIT and
            complete["rms_sample_max"] < COMMAND_RMS_LIMIT
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--controller", type=Path, required=True)
    parser.add_argument("--expect-source-sha256")
    parser.add_argument("--expect-controller-sha256")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check-limits", action="store_true")
    args = parser.parse_args()
    result = audit(args.source, args.controller, args.expect_source_sha256,
                   args.expect_controller_sha256)
    payload = json.dumps(result, indent=2) + "\n"
    if args.output is None:
        print(payload, end="")
    else:
        args.output.write_text(payload)
    if args.check_limits and not result["complete_controller_component_pass"]:
        raise SystemExit("The checked WaveStar controller replay gate failed")


if __name__ == "__main__":
    main()

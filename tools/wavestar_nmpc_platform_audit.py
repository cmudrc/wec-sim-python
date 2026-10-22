"""Replay identical WaveStar NMPC inputs to locate platform divergence.

Run this tool with two Python architectures against the same committed
fixtures, then compare the resulting NPZ arrays. Saved source pitch and
requests are diagnostic inputs here; this is not an autonomous parity gate.

From the repository root, run ``PYTHONPATH=. python
tools/wavestar_nmpc_platform_audit.py --output arm.npz`` in one environment.
In the other, add ``--output x86.npz --compare arm.npz``. The two environments
must use the same NumPy and SciPy versions.
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
from wecsim.dynamics.wavestar import WaveStarLinkage


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "tests/fixtures/matlab/wavestar_local_full_source/source.npz"
QP = ROOT / "tests/fixtures/matlab/wavestar_nmpc_independent_qp/source.npz"
MOMENT = (ROOT / "tests/fixtures/matlab/wavestar_nmpc_independent_predictor"
          / "python_moment.csv")


def replay(stroke, source_command, *, recursive):
    observer = WaveStarNmpcObserver()
    predictor = WaveStarNmpcPredictor()
    controller = WaveStarNmpcController()
    states = np.empty((len(stroke), 5))
    forecasts = np.empty((len(stroke), 40))
    commands = np.empty(len(stroke))
    previous = 0.
    for sample, position in enumerate(stroke):
        prior = previous if recursive else (source_command[sample - 1]
                                            if sample else 0.)
        state = observer.step(float(position), float(prior))
        forecast = predictor.step(float(state[-1]))
        command = controller.step(state, forecast, sample * .05)
        states[sample] = state
        forecasts[sample] = forecast
        commands[sample] = command
        previous = command
    return states, forecasts, commands


def audit(until):
    count = round(until / .05) + 1
    if not np.isclose((count - 1) * .05, until, rtol=0, atol=1e-10):
        raise ValueError("until must align with the 0.05 s controller clock")
    with np.load(SOURCE, allow_pickle=False) as source, \
            np.load(QP, allow_pickle=False) as qp:
        time = source["time"][:count].copy()
        pitch = source["pitch"][:count].copy()
        source_command = source["command"][:count].copy()
        qp_count = max(0, count - 300)
        qp_time = qp["time"][:qp_count].copy()
        qp_state = qp["estimated_state"][:qp_count].copy()
        qp_forecast = qp["forecast"][:, :qp_count].T.copy()
    if count > 4501 or until < 15 or not np.isfinite(pitch).all():
        raise ValueError("until must be between 15 and 225 s")
    np.testing.assert_allclose(time, np.arange(count) * .05, rtol=0, atol=1e-10)
    np.testing.assert_allclose(qp_time, time[300:], rtol=0, atol=1e-10)
    stroke = WaveStarLinkage().pto_stroke(pitch)
    fixed_state, fixed_forecast, fixed_command = replay(
        stroke, source_command, recursive=False,
    )
    recursive_state, recursive_forecast, recursive_command = replay(
        stroke, source_command, recursive=True,
    )
    qp_controller = WaveStarNmpcController()
    qp_command = np.array([
        qp_controller.step(state, forecast, float(instant))
        for state, forecast, instant in zip(qp_state, qp_forecast, qp_time)
    ])
    moment = np.loadtxt(MOMENT, delimiter=",")[:count]
    np.testing.assert_allclose(moment[:, 0], time, rtol=0, atol=1e-10)
    predictor = WaveStarNmpcPredictor()
    predictor_forecast = np.array([
        predictor.step(float(value)) for value in moment[:, 1]
    ])
    return {
        "architecture": np.array(platform.machine()),
        "numpy_version": np.array(np.__version__),
        "scipy_version": np.array(scipy.__version__),
        "fixture_hashes": np.array([
            sha256(path.read_bytes()).hexdigest()
            for path in (SOURCE, QP, MOMENT)
        ]),
        "time": time,
        "stroke": stroke,
        "source_command": source_command,
        "fixed_state": fixed_state,
        "fixed_forecast": fixed_forecast,
        "fixed_command": fixed_command,
        "recursive_state": recursive_state,
        "recursive_forecast": recursive_forecast,
        "recursive_command": recursive_command,
        "qp_command": qp_command,
        "predictor_forecast": predictor_forecast,
    }


def compare(first, second):
    """Summarize component gaps on identical diagnostic inputs."""
    for key in ("time", "source_command", "fixture_hashes"):
        np.testing.assert_array_equal(first[key], second[key])
    # The pinned pitch is identical, but linkage arithmetic can differ by a
    # few ulps across architectures before the controller runs.
    np.testing.assert_allclose(first["stroke"], second["stroke"],
                               rtol=0, atol=1e-15)
    for key in ("numpy_version", "scipy_version"):
        if first[key].item() != second[key].item():
            raise ValueError(f"Different {key} values cannot isolate architecture")
    time = first["time"]
    stroke_gap = np.abs(first["stroke"] - second["stroke"])
    summary = {"stroke_input": {
        "max_abs": float(stroke_gap.max()),
        "max_time_s": float(time[np.argmax(stroke_gap)]),
    }}
    for key in (
        "fixed_state", "fixed_forecast", "fixed_command",
        "recursive_state", "recursive_forecast", "recursive_command",
        "qp_command", "predictor_forecast",
    ):
        left = first[key]
        right = second[key]
        if left.shape != right.shape or not np.isfinite(left).all() \
                or not np.isfinite(right).all():
            raise ValueError(f"Invalid or mismatched {key} arrays")
        times = time[300:] if key == "qp_command" else time
        selected = times >= 15
        difference = np.abs(left[selected] - right[selected])
        by_sample = difference.reshape(len(difference), -1).max(axis=1)
        peak = int(np.argmax(by_sample))
        above = np.flatnonzero(by_sample > 1e-6)
        summary[key] = {
            "max_abs": float(by_sample[peak]),
            "max_time_s": float(times[selected][peak]),
            "first_above_1e-6_s": (float(times[selected][above[0]])
                                    if len(above) else None),
        }
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--until", type=float, default=75.)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--compare", type=Path,
                        help="compare this run with another architecture's NPZ")
    args = parser.parse_args()
    results = audit(args.until)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.output, **results)
    report = {
        "output": str(args.output),
        "architecture": platform.machine(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "samples": len(results["time"]),
    }
    if args.compare:
        with np.load(args.compare, allow_pickle=False) as other:
            report["comparison"] = compare(other, results)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()

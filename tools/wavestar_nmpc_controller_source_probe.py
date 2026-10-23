"""Replay logged MATLAB NMPC inputs to locate closed-loop divergence.

This is diagnostic only. Each Python controller sees the corresponding source
observer states and AR forecast, never an independently integrated plant.
"""

import argparse
import json
from hashlib import sha256
from pathlib import Path

import numpy as np
from scipy.io import loadmat

from wecsim.control.wavestarNmpc import WaveStarNmpcController


def _values(log, key, count):
    signal = log[key]
    time = np.asarray(signal["time"]).reshape(-1)
    np.testing.assert_allclose(time[:count], np.arange(count) * .05,
                               rtol=0, atol=1e-9)
    values = np.asarray(signal["signals"]["values"])
    if values.shape[-1] == len(time) and values.ndim == 2:
        values = values.T
    if values.shape[0] < count:
        raise ValueError(f"{key} has fewer than {count} samples")
    values = values[:count]
    if not np.isfinite(values).all():
        raise ValueError(f"{key} contains nonfinite samples")
    return values


def replay(path, count):
    log = loadmat(path, simplify_cells=True)
    states = _values(log, "estimated_states", count)
    forecast = _values(log, "AR_excM_pred", count)
    command = _values(log, "cmd_ptoM", count).reshape(-1)
    if states.shape != (count, 5) or forecast.shape != (count, 40):
        raise ValueError("The source observer or AR forecast has an unexpected shape")
    controller = WaveStarNmpcController()
    replayed = np.array([
        controller.step(states[index], forecast[index], index * .05)
        for index in range(count)
    ])
    return states, forecast, command, replayed


def _maximum_gap(left, right):
    difference = np.abs(left - right)
    if difference.ndim == 2:
        difference = difference.max(axis=1)
    index = int(difference.argmax())
    return {"value": float(difference[index]), "time_s": round(index * .05, 10)}


def probe(source_a, source_b, end_time):
    count = round(end_time / .05) + 1
    if count < 1 or abs((count - 1) * .05 - end_time) > 1e-9:
        raise ValueError("The end time must align with the 50 ms controller")
    a = replay(source_a, count)
    b = replay(source_b, count)
    names = ("estimated_states", "AR_excM_pred", "source_command",
             "python_replayed_command")
    return {
        "end_time_s": end_time,
        "source_a_sha256": sha256(source_a.read_bytes()).hexdigest(),
        "source_b_sha256": sha256(source_b.read_bytes()).hexdigest(),
        "python_qp_replay_vs_source_command": {
            "a": _maximum_gap(a[2], a[3]),
            "b": _maximum_gap(b[2], b[3]),
        },
        "source_a_vs_b": {
            name: _maximum_gap(left, right)
            for name, left, right in zip(names, a, b)
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_a", type=Path)
    parser.add_argument("source_b", type=Path)
    parser.add_argument("--end-time", type=float, default=40.)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = probe(args.source_a, args.source_b, args.end_time)
    payload = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(payload, end="")
    else:
        args.output.write_text(payload)


if __name__ == "__main__":
    main()

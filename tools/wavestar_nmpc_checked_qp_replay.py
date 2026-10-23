"""Prepare and compare a direct MATLAB QP replay on checked source inputs."""

import argparse
import json
from hashlib import sha256
from pathlib import Path

import numpy as np
from scipy.io import loadmat, savemat

from wecsim import WaveStarNmpcController


def _checked(path, expected_hash):
    digest = sha256(path.read_bytes()).hexdigest()
    if digest != expected_hash:
        raise ValueError(f"{path} differs from its checked SHA-256")
    return digest


def _inputs(controller_path, matrices_path, controller_hash, matrices_hash):
    _checked(controller_path, controller_hash)
    _checked(matrices_path, matrices_hash)
    with np.load(controller_path, allow_pickle=False) as log:
        time = log["time"].copy()
        states = log["estimated_state"].copy()
        forecasts = log["forecast"].copy()
        requests = log["command"].copy()
    with np.load(matrices_path, allow_pickle=False) as saved:
        matrices = {key: saved[key].copy() for key in ("A", "Bu", "Bw", "R")}
    if (time.shape != (4501,) or states.shape != (4501, 5) or
            forecasts.shape != (4501, 40) or requests.shape != time.shape):
        raise ValueError("Checked source controller arrays have unexpected shapes")
    np.testing.assert_allclose(time, np.arange(4501) * .05, rtol=0, atol=1e-9)
    active = slice(300, None)
    controller = WaveStarNmpcController()
    for key, expected in (
        ("A", controller._transition),
        ("Bu", controller._input[:, None]),
        ("Bw", controller._wave_input[:, None]),
        ("R", controller._penalty),
    ):
        np.testing.assert_array_equal(matrices[key], expected)
    arrays = {
        "time": time[active],
        "states": states[active].T,
        "forecasts": forecasts[active].T,
        "source_commands": requests[active],
        **matrices,
    }
    if not all(np.isfinite(value).all() for value in arrays.values()):
        raise ValueError("Checked source controller arrays contain nonfinite values")
    return arrays


def _stats(delta, time):
    absolute = np.abs(delta)
    maximum = int(np.argmax(absolute))
    return {
        "max_abs_Nm": float(absolute[maximum]),
        "time_of_max_s": float(time[maximum]),
        "rms_Nm": float(np.sqrt(np.mean(np.square(delta)))),
        "max_15_to_20s_Nm": float(np.max(absolute[time <= 20.])),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("prepare", "compare"))
    parser.add_argument("--controller", type=Path, required=True)
    parser.add_argument("--matrices", type=Path, required=True)
    parser.add_argument("--controller-sha256", required=True)
    parser.add_argument("--matrices-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--matlab-result", type=Path)
    args = parser.parse_args()
    arrays = _inputs(args.controller, args.matrices, args.controller_sha256,
                     args.matrices_sha256)
    if args.mode == "prepare":
        matlab_arrays = {**arrays,
                         "time": arrays["time"][:, None],
                         "source_commands": arrays["source_commands"][:, None]}
        savemat(args.output, matlab_arrays, do_compression=True)
        return
    if args.matlab_result is None:
        parser.error("compare requires --matlab-result")
    saved = loadmat(args.matlab_result)
    time = np.asarray(saved["time"]).ravel()
    commands = np.asarray(saved["commands"]).ravel()
    plans = np.asarray(saved["plans"])
    np.testing.assert_array_equal(time, arrays["time"])
    if commands.shape != time.shape or plans.shape != (len(time), 40):
        raise ValueError("MATLAB QP replay has unexpected shapes")
    if not np.isfinite(commands).all() or not np.isfinite(plans).all():
        raise ValueError("MATLAB QP replay contains nonfinite values")
    np.testing.assert_allclose(commands, plans[:, 0], rtol=0, atol=1e-12)

    controller = WaveStarNmpcController()
    python_plans = np.empty_like(plans)
    for index, instant in enumerate(time):
        controller.step(arrays["states"][:, index],
                        arrays["forecasts"][:, index], float(instant))
        python_plans[index] = controller._nominal
    result = {
        "samples": len(time),
        "source_controller_sha256": args.controller_sha256,
        "matrices_sha256": args.matrices_sha256,
        "matlab_replay_sha256": sha256(args.matlab_result.read_bytes()).hexdigest(),
        "source_vs_matlab_qp": _stats(arrays["source_commands"] - commands, time),
        "python_vs_matlab_qp": _stats(python_plans[:, 0] - commands, time),
        "python_vs_source": _stats(python_plans[:, 0] - arrays["source_commands"], time),
        "full_plan_python_vs_matlab_max_abs_Nm": float(
            np.max(np.abs(python_plans - plans))
        ),
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

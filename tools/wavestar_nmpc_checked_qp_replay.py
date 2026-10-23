"""Prepare and compare a direct MATLAB QP replay on checked source inputs."""

import argparse
import json
import platform
from hashlib import sha256
from pathlib import Path

import numpy as np
import scipy
from scipy.io import loadmat, savemat

from wecsim import WaveStarNmpcController


COMMAND_MAX_LIMIT = .01
COMMAND_RMS_LIMIT = .0003


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
        # SciPy's matrix exponential can differ by a few ulps across CPUs.
        # Both QP replays use the checked matrix, not the local recomputation.
        np.testing.assert_allclose(matrices[key], expected, rtol=0, atol=1e-14)
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
    parser.add_argument("--matlab-replay-sha256")
    parser.add_argument("--allow-source-difference", action="store_true")
    parser.add_argument("--check-limits", action="store_true")
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
    if args.matlab_replay_sha256 is not None:
        _checked(args.matlab_result, args.matlab_replay_sha256)
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
    source_replay_gap = _stats(arrays["source_commands"] - commands, time)
    if (not args.allow_source_difference and
            source_replay_gap["max_abs_Nm"] >= 1e-8):
        raise ValueError("Direct MATLAB QP did not reproduce the source requests")

    controller = WaveStarNmpcController()
    controller._transition = arrays["A"].copy()
    controller._input = arrays["Bu"][:, 0].copy()
    controller._wave_input = arrays["Bw"][:, 0].copy()
    frozen_prior = WaveStarNmpcController()
    frozen_prior._transition = arrays["A"].copy()
    frozen_prior._input = arrays["Bu"][:, 0].copy()
    frozen_prior._wave_input = arrays["Bw"][:, 0].copy()
    python_plans = np.empty_like(plans)
    frozen_prior_plans = np.empty_like(plans)
    for index, instant in enumerate(time):
        controller.step(arrays["states"][:, index],
                        arrays["forecasts"][:, index], float(instant))
        python_plans[index] = controller._nominal
        frozen_prior._nominal = (plans[index - 1].copy() if index else
                                 np.zeros(40))
        frozen_prior.step(arrays["states"][:, index],
                          arrays["forecasts"][:, index], float(instant))
        frozen_prior_plans[index] = frozen_prior._nominal
    one_step_first_gap = frozen_prior_plans[:, 0] - commands
    peak_one_step = int(np.argmax(np.abs(one_step_first_gap)))
    result = {
        "samples": len(time),
        "python_platform": platform.platform(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "source_controller_sha256": args.controller_sha256,
        "matrices_sha256": args.matrices_sha256,
        "matlab_replay_sha256": sha256(args.matlab_result.read_bytes()).hexdigest(),
        "source_vs_matlab_qp": source_replay_gap,
        "source_reproduction_limit_Nm": 1e-8,
        "source_reproduction_pass": source_replay_gap["max_abs_Nm"] < 1e-8,
        "python_vs_matlab_qp": _stats(python_plans[:, 0] - commands, time),
        "python_vs_source": _stats(python_plans[:, 0] - arrays["source_commands"], time),
        "full_plan_python_vs_matlab_max_abs_Nm": float(
            np.max(np.abs(python_plans - plans))
        ),
        "python_one_step_from_matlab_prior_vs_matlab_qp": _stats(
            frozen_prior_plans[:, 0] - commands, time
        ),
        "full_plan_one_step_from_matlab_prior_max_abs_Nm": float(
            np.max(np.abs(frozen_prior_plans - plans))
        ),
        "largest_one_step_request_gap": {
            "time_s": float(time[peak_one_step]),
            "matlab_Nm": float(commands[peak_one_step]),
            "python_Nm": float(frozen_prior_plans[peak_one_step, 0]),
            "matlab_distance_to_torque_limit_Nm": float(
                controller.torque_limit - abs(commands[peak_one_step])
            ),
        },
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    request_gap = result["python_vs_matlab_qp"]
    if (args.check_limits and
            (request_gap["max_abs_Nm"] >= COMMAND_MAX_LIMIT or
             request_gap["rms_Nm"] >= COMMAND_RMS_LIMIT)):
        raise SystemExit("The checked QP request component gate failed")


if __name__ == "__main__":
    main()

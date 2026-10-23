"""Prepare and compare one fixed-input replay of the pinned MATLAB AR predictor."""

import argparse
import json
from hashlib import sha256
from pathlib import Path

import numpy as np
from scipy.io import loadmat, savemat


ROOT = Path(__file__).resolve().parents[1]
CONTROLLER = (ROOT / "tests/fixtures/matlab/"
              "wavestar_nmpc_reference_convergence_225s/"
              "controller_0p25ms_compatible_amd7763.npz")
CONTROLLER_SHA256 = "06d55c1972b6e9130fbeeea52d29cc6691fddfb7408c687c950ce9522c2498e3"


def checked_input():
    digest = sha256(CONTROLLER.read_bytes()).hexdigest()
    if digest != CONTROLLER_SHA256:
        raise ValueError("The checked WaveStar controller fixture changed")
    with np.load(CONTROLLER, allow_pickle=False) as data:
        time = data["time"].copy()
        state = data["estimated_state"].copy()
        logged = data["forecast"].copy()
    if (time.shape != (4501,) or state.shape != (4501, 5)
            or logged.shape != (4501, 40)
            or not all(np.isfinite(array).all()
                       for array in (time, state, logged))):
        raise ValueError("The checked WaveStar controller arrays are invalid")
    np.testing.assert_allclose(time, np.arange(4501) * .05,
                               rtol=0, atol=1e-9)
    return time, state[:, 4], logged, digest


def compare(result_path, other_path=None):
    time, _, logged, digest = checked_input()
    replay = np.asarray(loadmat(result_path)["forecasts"])
    if replay.shape != logged.shape or not np.isfinite(replay).all():
        raise ValueError("The MATLAB WaveStar predictor result is invalid")

    def stats(difference):
        by_sample = np.max(np.abs(difference), axis=1)
        peak = int(np.argmax(by_sample))
        active = by_sample[300:]
        active_peak = int(np.argmax(active)) + 300
        return {
            "max_abs_Nm": float(by_sample[peak]),
            "time_of_max_s": float(time[peak]),
            "max_active_after_15s_Nm": float(by_sample[active_peak]),
            "active_time_of_max_s": float(time[active_peak]),
        }

    report = {
        "controller_fixture_sha256": digest,
        "samples": len(time),
        "horizon": replay.shape[1],
        "matlab_replay_vs_checked_source": stats(replay - logged),
    }
    if other_path is not None:
        other = np.asarray(loadmat(other_path)["forecasts"])
        if other.shape != replay.shape or not np.isfinite(other).all():
            raise ValueError("The other MATLAB predictor result is invalid")
        report["matlab_replay_vs_other_platform"] = stats(replay - other)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="mode", required=True)
    prepare = sub.add_parser("prepare")
    prepare.add_argument("--output", type=Path, required=True)
    comparison = sub.add_parser("compare")
    comparison.add_argument("--matlab-result", type=Path, required=True)
    comparison.add_argument("--other-matlab-result", type=Path)
    comparison.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.mode == "prepare":
        time, moment, _, _ = checked_input()
        savemat(args.output, {"time": time[:, None],
                              "moment": moment[:, None]})
        print(args.output)
    else:
        result = compare(args.matlab_result, args.other_matlab_result)
        args.output.write_text(json.dumps(result, indent=2) + "\n")
        print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

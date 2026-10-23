"""Gate the corrected continuous-heading source on two full-duration steps."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.io import loadmat


# Chosen before evaluating the half-step source run. These limits exceed the
# prior three-seed source/Python 0.01 s gaps but still require small full-run
# changes after halving the source step.
LIMITS = {
    "wave_m": 1e-11,
    "yaw_rad": 1e-4,
    "yaw_speed_rad_s": 2e-5,
    "pto_torque_nm": 2.4,
    "pto_work_percent": 0.02,
}


def _load(path, dt):
    data = loadmat(path, squeeze_me=True)
    arrays = {key: np.asarray(data[key], dtype=float).ravel()
              for key in ("time", "yaw", "yawRate", "ptoForce", "ptoPower")}
    arrays["wave"] = np.asarray(data["wave"], dtype=float)[:, 1]
    length = round(250 / dt) + 1
    assert all(len(value) == length for value in arrays.values()), path
    assert np.max(np.abs(arrays["time"] - np.arange(length) * dt)) < 1e-10, path
    assert all(np.isfinite(value).all() for value in arrays.values()), path
    work = -np.trapezoid(arrays["ptoPower"], arrays["time"])
    assert work > 0 and abs(work - float(data["work"])) < 1e-8, path
    return arrays, work


def audit(source_dir, seed):
    files = {dt: source_dir / f"seed-{seed}-dt-{dt:g}.mat"
             for dt in (0.01, 0.005)}
    coarse, coarse_work = _load(files[0.01], 0.01)
    fine, fine_work = _load(files[0.005], 0.005)
    differences = {
        key: float(np.max(np.abs(coarse[key] - fine[key][::2])))
        for key in ("time", "wave", "yaw", "yawRate", "ptoForce")
    }
    assert differences["time"] < 1e-10
    work_change = 100 * (fine_work / coarse_work - 1)
    within_limits = bool(
        differences["wave"] < LIMITS["wave_m"]
        and differences["yaw"] < LIMITS["yaw_rad"]
        and differences["yawRate"] < LIMITS["yaw_speed_rad_s"]
        and differences["ptoForce"] < LIMITS["pto_torque_nm"]
        and abs(work_change) < LIMITS["pto_work_percent"]
    )
    return {
        "seed": seed,
        "limits": LIMITS,
        "source_sha256": {f"{dt:g}": hashlib.sha256(path.read_bytes()).hexdigest()
                          for dt, path in files.items()},
        "pto_work_j": {"0.01": coarse_work, "0.005": fine_work},
        "max_common_grid_difference": differences,
        "pto_work_change_percent": work_change,
        "within_limits": within_limits,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", required=True, type=Path)
    parser.add_argument("--seed", required=True, type=int, choices=(1, 2, 3))
    parser.add_argument("--output", type=Path)
    parser.add_argument("--enforce", action="store_true")
    args = parser.parse_args()
    report = audit(args.source_dir, args.seed)
    contents = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(contents)
    else:
        print(contents, end="")
    if args.enforce and not report["within_limits"]:
        raise SystemExit("Continuous-heading source did not meet step limits")


if __name__ == "__main__":
    main()

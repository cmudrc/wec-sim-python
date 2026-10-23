"""Summarize independent held-yaw step refinements from pinned MATLAB runs.

The source files come from ``oswec_passive_yaw_step_refinement.m``. Optional
Python input must be a self-advanced physical result; this tool never feeds
source motion, forces, or heading events into that run.
"""

import argparse
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np
from scipy.io import loadmat


STEPS = (0.01, 0.005, 0.0025)
# The earlier prescribed-heading diagnostic gates yaw, speed, and work at
# 0.003 rad, 0.0002 rad/s, and 0.1%. The PTO is a 120,000 N m s/rad damper,
# so its speed limit implies 24 N m of torque. Require the same event count
# and event times within one published 0.01 s sample before using a source
# step as a physical reference.
LIMITS = {
    "yaw_rad": 0.003,
    "yaw_speed_rad_s": 0.0002,
    "pto_torque_nm": 24,
    "pto_work_percent": 0.1,
    "heading_event_s": 0.01,
    "wave_m": 1e-11,
}


def _heading_events(time, yaw, directions, incident_direction=10, threshold=1):
    """Replay the published coefficient-selection rule on accepted yaw states."""
    previous = None
    headings = np.empty(len(yaw))
    for index, angle in enumerate(yaw):
        relative = incident_direction - np.degrees(angle)
        distance = (relative - directions + 180) % 360 - 180
        nearest = int(np.argmin(np.abs(distance)))
        if abs(distance[nearest]) <= threshold:
            previous = relative - distance[nearest]
        elif previous is None or abs(relative - previous) > threshold:
            previous = relative
        headings[index] = previous
    changes = np.flatnonzero(np.abs(np.diff(headings)) > 1e-5) + 1
    return time[changes]


def _load_source(path, dt):
    data = loadmat(path, squeeze_me=True)
    run = {key: np.asarray(data[key], dtype=float).ravel()
           for key in ("time", "yaw", "yawRate", "ptoForce", "ptoPower")}
    run["wave"] = np.asarray(data["wave"], dtype=float)[:, 1]
    size = round(250 / dt) + 1
    assert all(len(value) == size for value in run.values()), path
    assert np.max(np.abs(run["time"] - np.arange(size) * dt)) < 1e-10, path
    assert all(np.isfinite(value).all() for value in run.values()), path
    work = -np.trapezoid(run["ptoPower"], run["time"])
    assert work > 0 and abs(work - float(data["work"])) < 1e-8, path
    return run, work


def _event_comparison(first, second):
    same_count = len(first) == len(second)
    return {
        "same_heading_event_count": same_count,
        "max_heading_event_time_difference_s": (
            float(np.max(np.abs(first - second))) if same_count and len(first)
            else 0.0 if same_count else None
        ),
    }


def _within_limits(metrics):
    event_gap = metrics["max_heading_event_time_difference_s"]
    return (metrics["max_wave_difference"] < LIMITS["wave_m"]
            and metrics["max_yaw_difference"] < LIMITS["yaw_rad"]
            and metrics["max_yawRate_difference"] < LIMITS["yaw_speed_rad_s"]
            and metrics["max_ptoForce_difference"] < LIMITS["pto_torque_nm"]
            and abs(metrics["pto_work_difference_percent"]) < LIMITS["pto_work_percent"]
            and metrics["same_heading_event_count"]
            and event_gap is not None
            and event_gap <= LIMITS["heading_event_s"])


def audit(source_dir, hydro, seed, physical=None):
    with h5py.File(hydro) as file:
        directions = np.asarray(file["simulation_parameters/wave_dir"]).ravel()
    runs = {}
    summary = {"seed": seed, "limits": LIMITS,
               "matlab": {}, "source_step_pairs": {}}
    for dt in STEPS:
        path = source_dir / f"seed-{seed}-dt-{dt:g}.mat"
        run, work = _load_source(path, dt)
        events = _heading_events(run["time"], run["yaw"], directions)
        runs[dt] = (run, work)
        summary["matlab"][f"{dt:g}"] = {
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "pto_work_j": work,
            "heading_events_s": events.tolist(),
        }
    for coarse_dt, fine_dt in zip(STEPS, STEPS[1:], strict=False):
        coarse, coarse_work = runs[coarse_dt]
        fine, fine_work = runs[fine_dt]
        pair = {}
        for key in ("time", "wave", "yaw", "yawRate", "ptoForce"):
            pair[f"max_{key}_difference"] = float(np.max(
                np.abs(coarse[key] - fine[key][::2])))
        assert pair["max_time_difference"] < 1e-10
        assert pair["max_wave_difference"] < 1e-11
        pair["pto_work_difference_percent"] = 100 * (fine_work / coarse_work - 1)
        pair.update(_event_comparison(
            np.asarray(summary["matlab"][f"{coarse_dt:g}"]["heading_events_s"]),
            np.asarray(summary["matlab"][f"{fine_dt:g}"]["heading_events_s"]),
        ))
        pair["within_limits"] = _within_limits(pair)
        summary["source_step_pairs"][f"{coarse_dt:g}/{fine_dt:g}"] = pair
    if physical is not None:
        with np.load(physical) as archive:
            result = {key: archive[key] for key in archive.files}
        fine, fine_work = runs[STEPS[-1]]
        assert np.max(np.abs(result["time"] - fine["time"])) < 1e-10
        assert np.max(np.abs(result["wave"] - fine["wave"])) < 1e-11
        physical_events = _heading_events(
            result["time"], result["yaw"], directions,
        )
        physical_work = np.trapezoid(result["power"], result["time"])
        comparison = {
            "sha256": hashlib.sha256(physical.read_bytes()).hexdigest(),
            "pto_work_j": physical_work,
            "heading_events_s": physical_events.tolist(),
            "max_wave_difference": float(np.max(np.abs(result["wave"] - fine["wave"]))),
            "max_yaw_difference": float(np.max(np.abs(result["yaw"] - fine["yaw"]))),
            "max_yawRate_difference": float(np.max(np.abs(result["speed"] - fine["yawRate"]))),
            "max_ptoForce_difference": float(np.max(np.abs(result["torque"] - fine["ptoForce"]))),
            "pto_work_difference_percent": 100 * (physical_work / fine_work - 1),
        }
        comparison.update(_event_comparison(
            np.asarray(summary["matlab"][f"{STEPS[-1]:g}"]["heading_events_s"]),
            physical_events,
        ))
        comparison["within_limits"] = _within_limits(comparison)
        summary["physical_python_at_0.0025"] = comparison
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", required=True, type=Path)
    parser.add_argument("--hydro", required=True, type=Path)
    parser.add_argument("--seed", required=True, type=int, choices=(1, 2, 3))
    parser.add_argument("--physical", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    summary = audit(args.source_dir, args.hydro, args.seed, args.physical)
    contents = json.dumps(summary, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(contents)
    else:
        print(contents, end="")


if __name__ == "__main__":
    main()

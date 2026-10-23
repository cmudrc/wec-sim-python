"""Compare full-duration pinned WaveStar NMPC source step and repeat runs."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.io import loadmat


ROOT = Path(__file__).resolve().parents[1]
PINNED = ROOT / "tests/fixtures/matlab/wavestar_nmpc_split_step_225s/source.npz"
PREFIX = "WECCCOMP_NMPC_SPLIT_STEP_FULL_"
CASE = "WECCCOMP_Nonlinear_Model_Predictive"
ARTIFACTS = {
    "1ms": ("wavestar-nmpc-source-reference-1ms", .001),
    "0p5ms": ("wavestar-nmpc-source-reference-0p5ms", .0005),
    "0p25ms": ("wavestar-nmpc-source-reference-0p25ms", .00025),
}
FIELDS = {
    "pitch_rad": "pitch",
    "pitch_speed_rad_s": "pitch_speed",
    "stroke_m": "stroke",
    "axial_force_N": "axial_force",
    "requested_torque_Nm": "command",
}


def _read_artifact(path, step):
    audit = json.loads((path / "wavestar-reference-audit.json").read_text())
    if audit["source"]["plant_step_s"] != step or audit["source"]["samples"] != 4501:
        raise ValueError(f"Unexpected source provenance in {path}")
    body = np.loadtxt(path / f"{PREFIX}{CASE}_body1.csv", delimiter=",", ndmin=2)
    pto = np.loadtxt(path / f"{PREFIX}{CASE}_pto1.csv", delimiter=",", ndmin=2)
    log = loadmat(path / f"{PREFIX}controller.mat", simplify_cells=True)
    command = np.asarray(log["cmd_ptoM"]["signals"]["values"]).reshape(-1)
    time = np.arange(4501) * .05
    if body.shape != (4501, 25) or pto.shape != (4501, 43) \
            or command.shape != time.shape:
        raise ValueError(f"Unexpected source trace shape in {path}")
    for actual in (body[:, 0], pto[:, 0],
                   np.asarray(log["cmd_ptoM"]["time"]).reshape(-1)):
        np.testing.assert_allclose(actual, time, rtol=0, atol=1e-9)
    trace = {
        "pitch": body[:, 5],
        "pitch_speed": body[:, 11],
        "stroke": pto[:, 3],
        "axial_force": pto[:, 33],
        "command": command,
    }
    if not all(np.isfinite(value).all() for value in trace.values()):
        raise ValueError(f"Nonfinite source trace in {path}")
    return trace, audit["pointwise_limits"]


def _compare(left, right, limits):
    result = {}
    passed = True
    for label, key in FIELDS.items():
        difference = np.abs(right[key] - left[key])
        indices = np.flatnonzero(difference >= limits[label])
        peak = float(np.max(difference))
        result[label] = {
            "max_abs": peak,
            "first_limit_crossing_s": (round(float(indices[0] * .05), 10)
                                       if indices.size else None),
        }
        passed &= peak < limits[label]
    result["pitch_rad_first_15p2s"] = float(np.max(np.abs(
        right["pitch"][:305] - left["pitch"][:305],
    )))
    result["within_pointwise_limits"] = bool(passed)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    traces = {}
    limits = None
    for label, (artifact, step) in ARTIFACTS.items():
        trace, current_limits = _read_artifact(args.directory / artifact, step)
        if limits is not None and current_limits != limits:
            raise ValueError("The source audits used different pointwise limits")
        traces[label] = trace
        limits = current_limits
    with np.load(PINNED, allow_pickle=False) as baseline:
        older_1ms = {key: baseline[key].copy() for key in FIELDS.values()}
    result = {
        "pointwise_limits": limits,
        "pinned_vs_fresh_1ms": _compare(older_1ms, traces["1ms"], limits),
        "fresh_1ms_vs_0p5ms": _compare(traces["1ms"], traces["0p5ms"], limits),
        "fresh_0p5ms_vs_0p25ms": _compare(
            traces["0p5ms"], traces["0p25ms"], limits,
        ),
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

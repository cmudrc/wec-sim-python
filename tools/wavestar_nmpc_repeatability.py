"""Compare two pinned MATLAB WaveStar NMPC runs sharing one hydro fit."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.io import loadmat


PREFIX = "WECCCOMP_NMPC_SOURCE_"
CASE = "WECCCOMP_Nonlinear_Model_Predictive_"


def csv(directory: Path, name: str) -> np.ndarray:
    values = np.loadtxt(directory / f"{PREFIX}{name}.csv", delimiter=",", ndmin=2)
    if not np.isfinite(values).all():
        raise ValueError(f"Nonfinite MATLAB values in {directory}/{name}")
    return values


def command(directory: Path) -> tuple[np.ndarray, np.ndarray]:
    saved = loadmat(directory / f"{PREFIX}controller.mat", simplify_cells=True)
    signal = saved["cmd_ptoM"]
    return np.asarray(signal["time"]).ravel(), np.asarray(signal["signals"]["values"]).ravel()


def difference(time: np.ndarray, first: np.ndarray, second: np.ndarray) -> dict:
    if first.shape != second.shape or first.shape != time.shape:
        raise ValueError(f"Source trace dimensions differ: {first.shape}, {second.shape}, {time.shape}")
    if not np.isfinite(first).all() or not np.isfinite(second).all():
        raise ValueError("Nonfinite MATLAB source trace")
    error = np.abs(first - second)
    peak = int(np.argmax(error))
    return {
        "max_abs": float(error[peak]),
        "max_time_s": float(time[peak]),
        "rms": float(np.sqrt(np.mean(error**2))),
        "final_abs": float(error[-1]),
        "first_above_1e-5_s": next((float(time[i]) for i in np.flatnonzero(error > 1e-5)), None),
    }


def compare(first_dir: Path, second_dir: Path) -> dict:
    first_wave = csv(first_dir, "wave")
    second_wave = csv(second_dir, "wave")
    first_components = csv(first_dir, "components")
    second_components = csv(second_dir, "components")
    if first_wave.shape != (4501, 2) or first_components.shape != (500, 4):
        raise ValueError("Unexpected published WaveStar NMPC sea dimensions")
    np.testing.assert_allclose(first_wave, second_wave, rtol=0, atol=1e-12)
    np.testing.assert_allclose(first_components, second_components, rtol=0, atol=1e-12)

    first_body = csv(first_dir, CASE + "body1")
    second_body = csv(second_dir, CASE + "body1")
    first_pto = csv(first_dir, CASE + "pto1")
    second_pto = csv(second_dir, CASE + "pto1")
    if first_body.shape != (4501, 25) or first_pto.shape != (4501, 43):
        raise ValueError("Unexpected published WaveStar NMPC body/PTO dimensions")
    for left, right in ((first_body, second_body), (first_pto, second_pto)):
        np.testing.assert_allclose(left[:, 0], right[:, 0], rtol=0, atol=1e-12)
    first_time, first_command = command(first_dir)
    second_time, second_command = command(second_dir)
    np.testing.assert_allclose(first_time, second_time, rtol=0, atol=1e-12)
    np.testing.assert_allclose(first_time, first_body[:, 0], rtol=0, atol=1e-10)

    time = first_body[:, 0]
    signals = {
        "float_pitch_rad": (first_body[:, 5], second_body[:, 5]),
        "float_pitch_speed_rad_s": (first_body[:, 11], second_body[:, 11]),
        "pto_stroke_m": (first_pto[:, 3], second_pto[:, 3]),
        "pto_actuation_N": (first_pto[:, 33], second_pto[:, 33]),
        "controller_command_Nm": (first_command, second_command),
    }
    windows = {"pre_control": time < 10, "resistive": (time >= 10) & (time < 15),
               "nmpc": time >= 15, "full": np.ones(time.shape, dtype=bool)}
    return {
        "source": "Pinned R2025b WEC-Sim and Applications, two MATLAB invocations, one generated wavestar.h5",
        "sea_identical_within": 1e-12,
        "windows": {
            name: {signal: difference(time[mask], left[mask], right[mask])
                   for signal, (left, right) in signals.items()}
            for name, mask in windows.items()
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("first", type=Path)
    parser.add_argument("second", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    result = compare(args.first, args.second)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

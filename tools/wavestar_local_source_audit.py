"""Compare a local unchanged WaveStar MATLAB run to the pinned hosted source."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.io import loadmat


ROOT = Path(__file__).resolve().parents[1]
BASELINE = ROOT / "tests/fixtures/matlab/wavestar_nmpc_full_step/baseline.npz"
HYDRO = ROOT / "tests/fixtures/matlab/wavestar_fixed_input_step/wavestar.h5"
HYDRO_SHA256 = "22702b94f052dfe6bb93d257006654d4e51246020a53f25b432d591bf8139d1a"


def _read_csv(folder: Path, name: str, shape: tuple[int, int]) -> np.ndarray:
    values = np.loadtxt(folder / f"{name}.csv", delimiter=",", ndmin=2)
    assert values.shape == shape and np.isfinite(values).all(), name
    return values


def _difference(time: np.ndarray, local: np.ndarray, hosted: np.ndarray) -> dict:
    gap = np.abs(local - hosted)
    return {
        "max_abs": float(gap.max()),
        "max_time_s": float(time[np.argmax(gap)]),
        "rms": float(np.sqrt(np.mean(gap**2))),
        "first_above_1e-4_s": next(
            (float(time[index]) for index in np.flatnonzero(gap > 1e-4)), None
        ),
    }


def audit(folder: Path) -> dict:
    assert hashlib.sha256(HYDRO.read_bytes()).hexdigest() == HYDRO_SHA256
    copied_hydro = folder / "WECCCOMP/hydroData/wavestar.h5"
    assert hashlib.sha256(copied_hydro.read_bytes()).hexdigest() == HYDRO_SHA256
    local = {
        "body": _read_csv(folder, "body", (4501, 25)),
        "pto": _read_csv(folder, "pto", (4501, 43)),
        "components": _read_csv(folder, "components", (500, 4)),
        "wave": _read_csv(folder, "wave", (4501, 2)),
    }
    controller = loadmat(folder / "controller.mat", simplify_cells=True)
    command_time = np.asarray(controller["cmd_ptoM"]["time"]).ravel()
    command = np.asarray(controller["cmd_ptoM"]["signals"]["values"]).ravel()
    with np.load(BASELINE, allow_pickle=False) as reference:
        time = reference["body"][:, 0]
        np.testing.assert_allclose(time, np.arange(4501) * .05, rtol=0, atol=1e-10)
        for name in ("body", "pto", "wave"):
            np.testing.assert_allclose(local[name][:, 0], time, rtol=0, atol=1e-10)
        np.testing.assert_allclose(command_time, time, rtol=0, atol=1e-10)
        assert command.shape == time.shape and np.isfinite(command).all()

        # These identity gates were fixed before inspecting the local run.
        np.testing.assert_allclose(local["components"], reference["components"],
                                   rtol=0, atol=1e-12)
        np.testing.assert_allclose(local["wave"], reference["wave"],
                                   rtol=0, atol=1e-12)
        early = time <= 20
        for name, columns, limit in (
            ("body", slice(1, 13), 1e-5),
            ("pto", slice(1, 13), 1e-5),
            ("pto", slice(31, 37), 1e-3),
        ):
            assert np.max(np.abs(local[name][early, columns]
                                 - reference[name][early, columns])) < limit
        assert np.max(np.abs(command[early] - reference["command"][early])) < 1e-4

        # Closed-loop trajectories can amplify platform-level differences.
        # Report all 225 s without treating this self-comparison as Python parity.
        signals = {
            "float_pitch_rad": (local["body"][:, 5], reference["body"][:, 5]),
            "float_pitch_speed_rad_s": (local["body"][:, 11], reference["body"][:, 11]),
            "pto_stroke_m": (local["pto"][:, 3], reference["pto"][:, 3]),
            "pto_actuation_N": (local["pto"][:, 33], reference["pto"][:, 33]),
            "controller_command_Nm": (command, reference["command"]),
        }
        windows = {"pre_control": time < 10,
                   "resistive": (time >= 10) & (time < 15),
                   "nmpc": time >= 15, "full": np.ones_like(time, dtype=bool)}
        return {
            "source": "Local pinned R2025b with installed Statistics/System Identification toolboxes",
            "input_identity_gate": "passed",
            "early_20s_gate": "passed",
            "windows": {
                label: {name: _difference(time[mask], lhs[mask], rhs[mask])
                        for name, (lhs, rhs) in signals.items()}
                for label, mask in windows.items()
            },
        }


def independent_python_audit(folder: Path) -> dict:
    """Compare one independently integrated Python trajectory with both sources."""
    from wecsim import WaveStarNmpcPTO
    from wecsim.dynamics.wavestar import run_wavestar_published
    from wecsim.environment.irregularWave import jonswap_equal_energy_components

    local_body = _read_csv(folder, "body", (4501, 25))
    local_pto = _read_csv(folder, "pto", (4501, 43))
    local_components = _read_csv(folder, "components", (500, 4))
    command_mat = loadmat(folder / "controller.mat", simplify_cells=True)
    local_command = np.asarray(
        command_mat["cmd_ptoM"]["signals"]["values"]
    ).ravel()
    sea = jonswap_equal_energy_components(
        HYDRO, significant_height=.1042, peak_period=1.836,
        directions=np.array([0.]), spreading=np.array([1.]), gamma=3.3,
        phase=local_components[:, 3, None],
    )
    controller = WaveStarNmpcPTO(plant_dt=.001, control_dt=.05)
    response = run_wavestar_published(
        HYDRO, sea, dt=.001, end_time=225, ramp_time=25,
        g=9.80665, pto_controller=controller, output_stride=50,
    )
    time = local_body[:, 0]
    np.testing.assert_allclose(response.time, time, rtol=0, atol=1e-9)
    measured = {
        "float_pitch_rad": response.angle,
        "float_pitch_speed_rad_s": response.angular_speed,
        "pto_stroke_m": response.pto_stroke,
        "pto_actuation_N": response.pto_force,
        "controller_command_Nm": controller.command_torque,
    }
    local = {
        "float_pitch_rad": local_body[:, 5],
        "float_pitch_speed_rad_s": local_body[:, 11],
        "pto_stroke_m": local_pto[:, 3],
        "pto_actuation_N": local_pto[:, 33],
        "controller_command_Nm": local_command,
    }
    with np.load(BASELINE, allow_pickle=False) as reference:
        hosted = {
            "float_pitch_rad": reference["body"][:, 5],
            "float_pitch_speed_rad_s": reference["body"][:, 11],
            "pto_stroke_m": reference["pto"][:, 3],
            "pto_actuation_N": reference["pto"][:, 33],
            "controller_command_Nm": reference["command"],
        }
        windows = {"pre_control": time < 10,
                   "resistive": (time >= 10) & (time < 15),
                   "nmpc": time >= 15,
                   "full": np.ones_like(time, dtype=bool)}
        differences = {
            label: {
                name: {
                    "python_vs_local": _difference(time[mask], values[mask], local[name][mask]),
                    "python_vs_hosted": _difference(time[mask], values[mask], hosted[name][mask]),
                }
                for name, values in measured.items()
            }
            for label, mask in windows.items()
        }
        work = lambda force, speed: float(np.trapezoid(force * speed, time))
        return {
            "integration": "independent public WaveStar plant and NMPC, 1 ms plant, 50 ms controller",
            "differences": differences,
            "pto_work_J": {
                "python": work(response.pto_force, response.pto_speed),
                "local_matlab": work(local_pto[:, 33], local_pto[:, 9]),
                "hosted_matlab": work(reference["pto"][:, 33], reference["pto"][:, 9]),
            },
        }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_dir", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--independent-python", action="store_true")
    args = parser.parse_args()
    result = audit(args.source_dir)
    if args.independent_python:
        result["independent_python"] = independent_python_audit(args.source_dir)
    report = json.dumps(result, indent=2) + "\n"
    if args.output is not None:
        args.output.write_text(report)
    print(report)


if __name__ == "__main__":
    main()

"""Compare WaveStar NMPC source motion with refined Python plant steps.

The MATLAB source files come from the pinned WECCCOMP_NMPC_SOURCE baseline.
This is a diagnostic of the published 0.05 s ode8 trajectory, not a parity
gate or a replacement for the independent closed-loop test.
"""

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat

from wecsim import WaveStarNmpcPTO
from wecsim.environment.irregularWave import jonswap_equal_energy_components
from wecsim.dynamics.wavestar import run_wavestar_published


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hydro", type=Path, required=True)
    parser.add_argument("--reference-dir", type=Path, required=True)
    parser.add_argument("--plot", type=Path, required=True)
    parser.add_argument("--duration", type=float, choices=(30.0, 225.0),
                        default=30.0)
    parser.add_argument("--quarter-ms", action="store_true",
                        help="also compare a 0.25 ms Python plant step")
    args = parser.parse_args()

    prefix = "WECCCOMP_NMPC_SOURCE_"
    reference = args.reference_dir
    count = round(args.duration / .05) + 1
    components = np.loadtxt(reference / f"{prefix}components.csv", delimiter=",")
    source_float = np.loadtxt(
        reference / f"{prefix}WECCCOMP_Nonlinear_Model_Predictive_body1.csv",
        delimiter=",")[:count]
    source_pto = np.loadtxt(
        reference / f"{prefix}WECCCOMP_Nonlinear_Model_Predictive_pto1.csv",
        delimiter=",")[:count]
    source_controller = loadmat(reference / f"{prefix}controller.mat",
                                simplify_cells=True)
    source_command = np.asarray(
        source_controller["cmd_ptoM"]["signals"]["values"]
    )[:count]
    if (source_float.shape != (count, 25)
            or source_pto.shape != (count, 43)
            or source_command.shape != (count,)):
        raise ValueError("published NMPC source does not cover the requested duration")
    sea = jonswap_equal_energy_components(
        args.hydro, significant_height=.1042, peak_period=1.836,
        directions=np.array([0.]), spreading=np.array([1.]),
        gamma=3.3, phase=components[:, 3, None],
    )
    results = {}
    steps = (.001, .0005, .00025) if args.quarter_ms else (.001, .0005)
    for dt in steps:
        pto = WaveStarNmpcPTO(plant_dt=dt, control_dt=.05)
        response = run_wavestar_published(
            args.hydro, sea, dt=dt, end_time=args.duration, ramp_time=25,
            g=9.80665, pto_controller=pto,
            output_stride=round(.05 / dt),
        )
        if (response.time.shape != (count,)
                or np.max(np.abs(response.time - source_float[:, 0])) > 1e-9):
            raise ValueError("source and Python output grids do not align")
        results[dt] = {
            "pitch": response.angle,
            "speed": response.angular_speed,
            "command": pto.command_torque,
            "force": response.pto_force,
            "pto_speed": response.pto_speed,
        }

    time = source_float[:, 0]
    windows = (("before_control", 0, 10),
               ("resistive", 10, 15),
               ("nmpc", 15, args.duration + .001))
    summary = {}
    for label, start, stop in windows:
        selected = (time >= start) & (time < stop)
        matlab = {"pitch": source_float[:, 5],
                  "speed": source_float[:, 11],
                  "command": source_command,
                  "force": source_pto[:, 33]}
        summary[label] = {}
        for name, target in matlab.items():
            summary[label][name] = {
                "python_step_difference": float(np.max(np.abs(
                    results[.001][name][selected]
                    - results[.0005][name][selected]))),
                "matlab_difference_1ms": float(np.max(np.abs(
                    results[.001][name][selected] - target[selected]))),
            }
            if args.quarter_ms:
                summary[label][name].update({
                    "python_step_difference_quarter_half_ms": float(np.max(np.abs(
                        results[.00025][name][selected]
                        - results[.0005][name][selected]))),
                    "matlab_difference_quarter_ms": float(np.max(np.abs(
                        results[.00025][name][selected] - target[selected]))),
                })

    if args.duration == 225:
        # The published feedback trajectory is not pointwise repeatable
        # across regenerated radiation fits. These whole-run quantities test
        # whether its overall motion and energy exchange are closer than the
        # pointwise traces, without treating a diagnostic as a passing gate.
        rms = lambda values: float(np.sqrt(np.mean(np.square(values))))
        source_aggregate = {
            "pitch_rms": rms(source_float[:, 5]),
            "pitch_peak": float(np.max(np.abs(source_float[:, 5]))),
            "pitch_speed_rms": rms(source_float[:, 11]),
            "axial_force_rms": rms(source_pto[:, 33]),
            "axial_force_peak": float(np.max(np.abs(source_pto[:, 33]))),
            "requested_torque_rms": rms(source_command),
            "pto_work_J": float(np.trapezoid(
                source_pto[:, 33] * source_pto[:, 9], time,
            )),
        }
        summary["published_aggregate"] = {"source": source_aggregate}
        for dt in steps:
            trial = results[dt]
            measured = {
                "pitch_rms": rms(trial["pitch"]),
                "pitch_peak": float(np.max(np.abs(trial["pitch"]))),
                "pitch_speed_rms": rms(trial["speed"]),
                "axial_force_rms": rms(trial["force"]),
                "axial_force_peak": float(np.max(np.abs(trial["force"]))),
                "requested_torque_rms": rms(trial["command"]),
                "pto_work_J": float(np.trapezoid(
                    trial["force"] * trial["pto_speed"], time,
                )),
            }
            summary["published_aggregate"][f"python_{dt:g}"] = {
                "values": measured,
                "relative_error": {
                    name: abs(value - source_aggregate[name])
                    / abs(source_aggregate[name])
                    for name, value in measured.items()
                },
            }

    fig, axes = plt.subplots(3, 1, figsize=(9, 9), sharex=True,
                             constrained_layout=True)
    axes[0].plot(time, source_float[:, 5], color="#202632", lw=1.8,
                 label="MATLAB ode8, 50 ms")
    axes[0].plot(time, results[.001]["pitch"], color="#147d9d", lw=1.25,
                 label="Python plant, 1 ms")
    axes[0].plot(time, results[.0005]["pitch"], color="#db762e", lw=1,
                 ls="--", label="Python plant, 0.5 ms")
    if args.quarter_ms:
        axes[0].plot(time, results[.00025]["pitch"], color="#80568d", lw=1,
                     ls=":", label="Python plant, 0.25 ms")
    axes[0].set_ylabel("Float pitch (rad)")
    axes[0].set_title("WaveStar NMPC: published source versus independent Python",
                      loc="left", weight="bold")
    axes[0].legend(loc="upper left", fontsize=9, ncol=2, frameon=False)
    plot_steps = ((.001, "#147d9d", "-"), (.0005, "#db762e", "--"))
    if args.quarter_ms:
        plot_steps += ((.00025, "#80568d", ":"),)
    for dt, color, style in plot_steps:
        axes[1].semilogy(time, np.maximum(np.abs(
            results[dt]["pitch"] - source_float[:, 5]), 1e-9),
            color=color, ls=style, label=f"{dt * 1000:g} ms vs MATLAB")
        axes[2].plot(time, np.abs(results[dt]["command"] - source_command),
                     color=color, ls=style,
                     label=f"{dt * 1000:g} ms vs MATLAB")
    axes[1].semilogy(time, np.maximum(np.abs(
        results[.001]["pitch"] - results[.0005]["pitch"]), 1e-9),
        color="#689069", label="Python step difference")
    axes[2].plot(time, np.abs(
        results[.001]["command"] - results[.0005]["command"]),
        color="#689069", label="Python step difference")
    if args.quarter_ms:
        for axis, name in ((axes[1], "pitch"), (axes[2], "command")):
            axis.plot(time, np.abs(
                results[.0005][name] - results[.00025][name]),
                color="#80568d", ls=":", label="0.25 vs 0.5 ms")
    axes[1].set_ylabel("Absolute pitch difference (rad)")
    axes[2].set_ylabel("Torque-command difference (N m)")
    axes[2].set_xlabel("Time (s)")
    for axis in axes[1:]:
        axis.legend(loc="upper left", fontsize=9, ncol=3, frameon=False)
    for axis in axes:
        for instant in (10, 15):
            axis.axvline(instant, color="#8b8d92", lw=.8, ls=":")
        axis.grid(alpha=.18)
        axis.set_xlim(0, args.duration)
    args.plot.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.plot, dpi=170)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

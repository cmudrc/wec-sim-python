"""Audit a freshly generated, pinned WaveStar split-step MATLAB reference.

The MATLAB checkout supplies only the hydrodynamic file and wave phases to the
autonomous Python run. Saved source motion, force, stroke, and commands are
loaded after that run for comparison, never as feedback inputs.
"""

import argparse
import json
import sys
from hashlib import sha256
from pathlib import Path

import numpy as np
import scipy
from scipy.io import loadmat

from wecsim.dynamics.wavestar import run_wavestar_published
from wecsim.environment.irregularWave import jonswap_equal_energy_components


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tests"))
from test_wavestar_split_step_portable import _DelayedSampledPTO  # noqa: E402


FIXTURE = ROOT / "tests/fixtures/matlab/wavestar_nmpc_split_step_225s"
PREFIX = "WECCCOMP_NMPC_SPLIT_STEP_FULL_"
CASE = "WECCCOMP_Nonlinear_Model_Predictive"
HYDRO_SHA256 = "b2541ff01c319987d9157d1be9fbc9ec39a7660d0572d8e316c428b8aa3359d1"
POINTWISE_LIMITS = {
    "pitch_rad": 0.0005,
    "pitch_speed_rad_s": 0.005,
    "stroke_m": 0.0001,
    "axial_force_N": 1.0,
    "requested_torque_Nm": 0.1,
}
AGGREGATE_LIMIT = 0.005


def _trace(log, name):
    signal = log[name]
    time = np.asarray(signal["time"]).reshape(-1)
    values = np.asarray(signal["signals"]["values"]).reshape(len(time), -1)
    return time, values


def _rms(values):
    return float(np.sqrt(np.mean(np.square(values))))


def audit(source_dir, hydro, step):
    digest = sha256(hydro.read_bytes()).hexdigest()
    if digest != HYDRO_SHA256:
        raise ValueError(f"HDF5 hash differs from the pinned fixture: {digest}")
    with np.load(FIXTURE / "source.npz", allow_pickle=False) as fixture:
        pinned = {name: fixture[name].copy() for name in fixture.files}

    source_step = np.loadtxt(source_dir / f"{PREFIX}steps.csv", delimiter=",")
    np.testing.assert_allclose(source_step, [step, .05, .05], rtol=0, atol=1e-12)
    body = np.loadtxt(source_dir / f"{PREFIX}{CASE}_body1.csv",
                      delimiter=",", ndmin=2)
    pto = np.loadtxt(source_dir / f"{PREFIX}{CASE}_pto1.csv",
                     delimiter=",", ndmin=2)
    components = np.loadtxt(source_dir / f"{PREFIX}components.csv",
                            delimiter=",", ndmin=2)
    log = loadmat(source_dir / f"{PREFIX}controller.mat", simplify_cells=True)
    command_time, command_values = _trace(log, "cmd_ptoM")
    request = command_values[:, 0]
    time = np.arange(4501) * .05
    if body.shape != (4501, 25) or pto.shape != (4501, 43) \
            or components.shape != (500, 4) or request.shape != (4501,):
        raise ValueError("The full WaveStar source trace has unexpected shapes")
    for actual in (body[:, 0], pto[:, 0], command_time):
        np.testing.assert_allclose(actual, time, rtol=0, atol=1e-9)
    np.testing.assert_allclose(components, pinned["components"],
                               rtol=0, atol=1e-12)
    if not all(np.isfinite(array).all() for array in
               (body, pto, components, request)):
        raise ValueError("The source trace contains nonfinite samples")

    sea = jonswap_equal_energy_components(
        hydro, significant_height=.1042, peak_period=1.836,
        directions=np.array([0.]), spreading=np.array([1.]),
        gamma=3.3, phase=components[:, 3, None],
    )
    controller = _DelayedSampledPTO()
    response = run_wavestar_published(
        hydro, sea, dt=.001, end_time=225., ramp_time=25.,
        g=9.80665, pto_controller=controller, output_stride=50,
    )
    np.testing.assert_allclose(response.time, time, rtol=0, atol=1e-9)
    own_request = np.asarray(controller.commands)
    if own_request.shape != time.shape:
        raise ValueError("The Python controller omitted a sample")
    if not all(np.isfinite(array).all() for array in (
        response.angle, response.angular_speed, response.pto_stroke,
        response.pto_speed, response.pto_force, own_request,
    )):
        raise ValueError("The Python run contains nonfinite samples")

    pointwise = {
        "pitch_rad": float(np.max(np.abs(response.angle - body[:, 5]))),
        "pitch_speed_rad_s": float(np.max(np.abs(
            response.angular_speed - body[:, 11]))),
        "stroke_m": float(np.max(np.abs(response.pto_stroke - pto[:, 3]))),
        "axial_force_N": float(np.max(np.abs(response.pto_force - pto[:, 33]))),
        "requested_torque_Nm": float(np.max(np.abs(own_request - request))),
    }
    source_aggregate = {
        "pitch_rms": _rms(body[:, 5]),
        "pitch_peak": float(np.max(np.abs(body[:, 5]))),
        "pitch_speed_rms": _rms(body[:, 11]),
        "axial_force_rms": _rms(pto[:, 33]),
        "axial_force_peak": float(np.max(np.abs(pto[:, 33]))),
        "requested_torque_rms": _rms(request),
        "pto_work_J": float(np.trapezoid(pto[:, 33] * pto[:, 9], time)),
    }
    python_aggregate = {
        "pitch_rms": _rms(response.angle),
        "pitch_peak": float(np.max(np.abs(response.angle))),
        "pitch_speed_rms": _rms(response.angular_speed),
        "axial_force_rms": _rms(response.pto_force),
        "axial_force_peak": float(np.max(np.abs(response.pto_force))),
        "requested_torque_rms": _rms(own_request),
        "pto_work_J": float(np.trapezoid(
            response.pto_force * response.pto_speed, time,
        )),
    }
    if not all(np.isfinite(value) and abs(value) > 1e-10
               for value in source_aggregate.values()):
        raise ValueError("The source aggregates cannot form relative gates")
    aggregate_gap = {
        name: abs(python_aggregate[name] - value) / abs(value)
        for name, value in source_aggregate.items()
    }
    prefix_count = round(15.2 / .05) + 1
    baseline_gap = {
        "pitch_rad_first_15p2s": float(np.max(np.abs(
            body[:prefix_count, 5] - pinned["pitch"][:prefix_count]))),
        "pitch_rad_full": float(np.max(np.abs(body[:, 5] - pinned["pitch"]))),
        "requested_torque_Nm_full": float(np.max(np.abs(
            request - pinned["command"]))),
    }
    return {
        "source": {
            "wec_sim_commit": "0753b2e47f2457c078751dcfe5d251d1767b80ab",
            "applications_commit": "d53d4d4c9eda2581f04204f5d394a6ef84bb099e",
            "hydro_sha256": digest,
            "plant_step_s": step,
            "controller_step_s": .05,
            "samples": len(time),
        },
        "python": {"numpy": np.__version__, "scipy": scipy.__version__},
        "source_vs_pinned_1ms": baseline_gap,
        "python_physical_pointwise": pointwise,
        "pointwise_limits": POINTWISE_LIMITS,
        "pointwise_pass": all(pointwise[name] < limit
                              for name, limit in POINTWISE_LIMITS.items()),
        "python_physical_aggregate_relative_gap": aggregate_gap,
        "aggregate_limit": AGGREGATE_LIMIT,
        "aggregate_pass": all(value < AGGREGATE_LIMIT
                              for value in aggregate_gap.values()),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--hydro", type=Path, required=True)
    parser.add_argument("--step", type=float, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = audit(args.source_dir, args.hydro, args.step)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()

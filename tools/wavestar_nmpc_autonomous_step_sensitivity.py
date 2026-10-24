"""Audit 225 s WaveStar closed-loop sensitivity without source-driven feedback.

This is a negative diagnostic. It measures two pinned MATLAB source plant
steps and two independent Python physical plant steps against the unchanged
pointwise limits. It never supplies MATLAB states, strokes, forces, or torque
requests to the Python controller.
"""

import argparse
import json
import platform
import sys
from hashlib import sha256
from pathlib import Path

import numpy as np
import scipy

from wecsim.dynamics.wavestar import run_wavestar_published
from wecsim.environment.irregularWave import jonswap_equal_energy_components


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tests"))
from test_wavestar_split_step_portable import _DelayedSampledPTO  # noqa: E402


SOURCE = ROOT / "tests/fixtures/matlab/wavestar_nmpc_reference_convergence_225s"
HYDRO = ROOT / "tests/fixtures/matlab/wavestar_nmpc_split_step_225s/wavestar.h5"
FILES = {
    "source_0p125ms_macos_arm_m1.npz":
        "ee791a11bb3c199e510dcc33e9f1b718313a34473694258ce5f75fd44ae55df0",
    "source_0p0625ms_macos_arm_m1.npz":
        "b3fe52fe66256a218c924c1c039f60c4b3d9d2192f36b1d9d4efe56831fac26b",
    "wavestar.h5":
        "b2541ff01c319987d9157d1be9fbc9ec39a7660d0572d8e316c428b8aa3359d1",
}
LIMITS = {
    "pitch_rad": .0005,
    "pitch_speed_rad_s": .005,
    "stroke_m": .0001,
    "axial_force_N": 1.,
    "requested_torque_Nm": .1,
}
CHANNELS = {
    "pitch_rad": "pitch",
    "pitch_speed_rad_s": "pitch_speed",
    "stroke_m": "stroke",
    "axial_force_N": "axial_force",
    "requested_torque_Nm": "command",
}


def _checked(path, name):
    actual = sha256(path.read_bytes()).hexdigest()
    if actual != FILES[name]:
        raise ValueError(f"Pinned WaveStar input {path} changed: {actual}")
    return actual


def _load_source(name):
    path = SOURCE / name
    _checked(path, name)
    with np.load(path, allow_pickle=False) as saved:
        trace = {key: saved[key].copy() for key in saved.files}
    if trace["components"].shape != (500, 4):
        raise ValueError("The source wave components changed shape")
    time = trace["time"]
    np.testing.assert_allclose(time, np.arange(4501) * .05,
                               rtol=0, atol=1e-9)
    for key in (*CHANNELS.values(), "pto_speed"):
        if trace[key].shape != time.shape or not np.isfinite(trace[key]).all():
            raise ValueError(f"The source {key} trace is invalid")
    return trace


def _compare(left, right):
    time = left["time"]
    np.testing.assert_array_equal(time, right["time"])
    pointwise = {}
    for label, key in CHANNELS.items():
        difference = np.abs(left[key] - right[key])
        first = np.flatnonzero(difference >= LIMITS[label])
        pointwise[label] = {
            "max_abs": float(np.max(difference)),
            "first_limit_crossing_s": float(time[first[0]]) if len(first) else None,
        }
    return {
        "pointwise": pointwise,
        "within_original_pointwise_limits": all(
            pointwise[label]["max_abs"] < limit for label, limit in LIMITS.items()
        ),
    }


def _independent_python(hydro, sea, plant_step):
    controller = _DelayedSampledPTO(plant_dt=plant_step)
    response = run_wavestar_published(
        hydro, sea, dt=plant_step, end_time=225., ramp_time=25.,
        g=9.80665, pto_controller=controller,
        output_stride=round(.05 / plant_step),
    )
    trace = {
        "time": response.time.copy(),
        "pitch": response.angle.copy(),
        "pitch_speed": response.angular_speed.copy(),
        "stroke": response.pto_stroke.copy(),
        "axial_force": response.pto_force.copy(),
        "command": np.asarray(controller.commands),
    }
    if any(value.shape != (4501,) or not np.isfinite(value).all()
           for value in trace.values()):
        raise ValueError("The independent Python trace is incomplete")
    return trace


def audit():
    _checked(HYDRO, "wavestar.h5")
    half = _load_source("source_0p125ms_macos_arm_m1.npz")
    quarter = _load_source("source_0p0625ms_macos_arm_m1.npz")
    np.testing.assert_array_equal(half["components"], quarter["components"])
    sea = jonswap_equal_energy_components(
        HYDRO, significant_height=.1042, peak_period=1.836,
        directions=np.array([0.]), spreading=np.array([1.]), gamma=3.3,
        phase=half["components"][:, 3, None],
    )
    python_1ms = _independent_python(HYDRO, sea, .001)
    python_0p5ms = _independent_python(HYDRO, sea, .0005)
    return {
        "purpose": "negative full-duration autonomous pointwise diagnostic",
        "duration_s": 225.,
        "controller_step_s": .05,
        "python_feedback_inputs": "pinned HDF5, wave phases, and own physical states only",
        "original_pointwise_limits": LIMITS,
        "pinned_sha256": FILES,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "matlab_0p125ms_vs_0p0625ms": _compare(half, quarter),
        "python_1ms_vs_0p5ms": _compare(python_1ms, python_0p5ms),
        "python_1ms_vs_matlab_0p125ms": _compare(python_1ms, half),
        "python_0p5ms_vs_matlab_0p125ms": _compare(python_0p5ms, half),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = audit()
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

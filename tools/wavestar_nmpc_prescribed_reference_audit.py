"""Audit the physical WaveStar plant on delayed source torque requests.

The source request is the only motion-related runtime input. This component
diagnostic does not count as an independent autonomous closed-loop gate.
"""

import argparse
import json
import platform
from hashlib import sha256
from pathlib import Path

import numpy as np
import scipy

from wecsim import WaveStarNmpcActuator
from wecsim.dynamics.wavestar import run_wavestar_published
from wecsim.environment.irregularWave import jonswap_equal_energy_components


POINTWISE_LIMITS = {
    "pitch_rad": .0005,
    "pitch_speed_rad_s": .005,
    "stroke_m": .0001,
    "axial_force_N": 1.,
}


class _DelayedSourceTorque:
    """Apply the prior 50 ms request, matching the split-step source delay."""

    def __init__(self, command, plant_dt):
        self.plant_dt = plant_dt
        self.samples_per_tick = round(.05 / plant_dt)
        if (self.samples_per_tick < 1 or
                abs(self.samples_per_tick * plant_dt - .05) > 1e-12):
            raise ValueError("The plant step must divide 50 ms")
        self.command = command
        self.actuator = WaveStarNmpcActuator()
        self.sample = 0

    def step(self, stroke, *, noise=0., dropout=False):
        if noise != 0 or dropout:
            raise ValueError("This diagnostic uses noiseless source requests")
        if self.sample % self.samples_per_tick == 0:
            tick = self.sample // self.samples_per_tick
            prior = self.command[tick - 1] if tick else 0.
            force = self.actuator.step(prior, stroke)
        else:
            force = self.actuator.force_at(stroke)
        self.sample += 1
        return force

    def force_at_stroke(self, stroke):
        return self.actuator.force_at(stroke)


def audit(source_path, hydro_path, plant_step, end_time,
          expected_source_sha256=None, expected_hydro_sha256=None):
    samples = round(end_time / .05) + 1
    if samples < 2 or abs((samples - 1) * .05 - end_time) > 1e-9:
        raise ValueError("The end time must align with 50 ms")
    source_digest = sha256(source_path.read_bytes()).hexdigest()
    hydro_digest = sha256(hydro_path.read_bytes()).hexdigest()
    if (expected_source_sha256 is not None and
            source_digest != expected_source_sha256):
        raise ValueError("The source trace hash differs from the checked input")
    if (expected_hydro_sha256 is not None and
            hydro_digest != expected_hydro_sha256):
        raise ValueError("The hydrodynamic hash differs from the checked input")
    with np.load(source_path, allow_pickle=False) as source:
        saved = {key: source[key].copy() for key in source.files}
    required = ("time", "components", "pitch", "pitch_speed", "stroke",
                "pto_speed", "axial_force", "command")
    if any(key not in saved for key in required):
        raise ValueError("The source trace is incomplete")
    if any(saved[key].shape != (4501,) for key in required[2:]):
        raise ValueError("The source trace needs 4,501 samples")
    if saved["time"].shape != (4501,) or saved["components"].shape != (500, 4):
        raise ValueError("The source time or sea has an unexpected shape")
    if not all(np.isfinite(saved[key]).all() for key in required):
        raise ValueError("The source contains nonfinite samples")
    np.testing.assert_allclose(saved["time"][:samples],
                               np.arange(samples) * .05, rtol=0, atol=1e-9)

    sea = jonswap_equal_energy_components(
        hydro_path, significant_height=.1042, peak_period=1.836,
        directions=np.array([0.]), spreading=np.array([1.]), gamma=3.3,
        phase=saved["components"][:, 3, None],
    )
    stride = round(.05 / plant_step)
    controller = _DelayedSourceTorque(saved["command"], plant_step)
    response = run_wavestar_published(
        hydro_path, sea, dt=plant_step, end_time=end_time,
        ramp_time=25., g=9.80665, pto_controller=controller,
        output_stride=stride,
    )
    np.testing.assert_allclose(response.time, saved["time"][:samples],
                               rtol=0, atol=1e-9)
    fields = {
        "pitch_rad": (response.angle, saved["pitch"][:samples]),
        "pitch_speed_rad_s": (response.angular_speed,
                              saved["pitch_speed"][:samples]),
        "stroke_m": (response.pto_stroke, saved["stroke"][:samples]),
        "axial_force_N": (response.pto_force,
                          saved["axial_force"][:samples]),
    }
    pointwise = {}
    for name, (own, reference) in fields.items():
        difference = np.abs(own - reference)
        crossing = np.flatnonzero(difference >= POINTWISE_LIMITS[name])
        pointwise[name] = {
            "max_abs": float(difference.max()),
            "first_limit_crossing_s": (float(response.time[crossing[0]])
                                       if crossing.size else None),
        }
    source_work = float(np.trapezoid(
        saved["axial_force"][:samples] * saved["pto_speed"][:samples],
        saved["time"][:samples],
    ))
    python_work = float(np.trapezoid(
        response.pto_force * response.pto_speed, response.time,
    ))
    if abs(source_work) < 1e-10:
        raise ValueError("Source work is too small for a relative audit")
    return {
        "source_trace_sha256": source_digest,
        "hydro_sha256": hydro_digest,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "physical_plant_step_s": plant_step,
        "end_time_s": end_time,
        "prescribed_source_torque_with_one_sample_delay": True,
        "pointwise": pointwise,
        "pointwise_limits": POINTWISE_LIMITS,
        "pointwise_pass": all(
            pointwise[name]["max_abs"] < limit
            for name, limit in POINTWISE_LIMITS.items()
        ),
        "source_work_J": source_work,
        "python_work_J": python_work,
        "relative_work_gap": abs(python_work - source_work) / abs(source_work),
        "work_relative_limit": .005,
        "work_pass": abs(python_work - source_work) / abs(source_work) < .005,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--hydro", type=Path, required=True)
    parser.add_argument("--step", type=float, default=.00025)
    parser.add_argument("--end-time", type=float, default=225.)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check-limits", action="store_true")
    parser.add_argument("--expect-source-sha256")
    parser.add_argument("--expect-hydro-sha256")
    args = parser.parse_args()
    result = audit(args.source, args.hydro, args.step, args.end_time,
                   args.expect_source_sha256, args.expect_hydro_sha256)
    payload = json.dumps(result, indent=2) + "\n"
    if args.output is None:
        print(payload, end="")
    else:
        args.output.write_text(payload)
    if args.check_limits and (not result["pointwise_pass"] or
                              not result["work_pass"]):
        raise SystemExit("The prescribed-command physical component gate failed")


if __name__ == "__main__":
    main()

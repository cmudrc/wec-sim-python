"""Run the published irregular OSWEC yaw input at a chosen solver step.

This is a local refinement diagnostic for issue #621, not a parity gate.
The source-equivalent ode4 path has been checked against MATLAB through
160 s at 0.00125/0.0025 s and through 250 s at 0.005 s. Its longer refined
traces remain surrogate predictions until MATLAB runs validate them.
"""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from wecsim import PMWave, WEC, WorldPoint


HYDRO = (Path(__file__).resolve().parents[1] / "tests/fixtures/matlab/"
         "oswec_multi_wave/oswec.h5")
HEADING_BANK = np.arange(-30, 30.25, 0.25)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--solver", choices=("implicit", "wecsim_ode4"),
                        required=True)
    parser.add_argument("--dt", type=float, required=True)
    parser.add_argument("--duration", type=float, default=600)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if (not np.isfinite(args.dt) or args.dt <= 0
            or not np.isfinite(args.duration) or args.duration <= 0):
        parser.error("dt and duration must be positive and finite")

    wec = WEC("OSWEC irregular selected-heading refinement")
    flap = wec.body("flap", HYDRO, mass=12_700, inertia=(1.85e6,) * 3,
                    passive_yaw=True, yaw_heading_bank=HEADING_BANK)
    wec.body("base", HYDRO, mass=999, inertia=(999,) * 3)
    yaw = wec.coordinate("yaw", flap.move(
        "yaw", pivot=WorldPoint(0, 0, -8.9),
    ))
    wec.rotational_pto("hinge", yaw, damping=120_000)
    result = wec.run(
        PMWave(2.5, 8, direction=10, seed=1, phase_generator="matlab"),
        dt=args.dt, end_time=args.duration, ramp_time=100,
        radiation_memory=40, solver=args.solver,
    )
    angle = result.bodies["flap"].position[:, 5]
    relative = 10 - np.degrees(angle)
    right = np.clip(np.searchsorted(HEADING_BANK, relative),
                    0, len(HEADING_BANK) - 1)
    left = np.maximum(right - 1, 0)
    # A tie selects the lower bank entry, as the source's first min index does.
    nearest = np.where(np.abs(HEADING_BANK[right] - relative)
                       < np.abs(HEADING_BANK[left] - relative), right, left)
    heading = HEADING_BANK[nearest]
    pto = result.ptos["hinge"]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.output, time=result.time, yaw=angle,
        speed=result.bodies["flap"].velocity[:, 5],
        torque=pto.force, power=pto.absorbed_power,
        wave=result.wave_elevation, heading=heading,
    )
    print(json.dumps({
        "output": str(args.output), "solver": args.solver,
        "dt_s": args.dt, "duration_s": args.duration,
        "samples": len(result.time),
        "heading_switches": int(np.count_nonzero(np.diff(heading))),
        "pto_work_J": float(np.trapezoid(pto.absorbed_power, result.time)),
        "hydro_sha256": hashlib.sha256(HYDRO.read_bytes()).hexdigest(),
    }))


if __name__ == "__main__":
    main()

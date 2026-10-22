"""Advance a published OSWEC held-yaw sea through the public physical solver.

Only the pinned HDF5 and MATLAB-compatible generated wave phases enter the
run. Source trajectories, forces, and heading events are never inputs.
"""

import argparse
from pathlib import Path

import numpy as np

from wecsim import PMWave, WEC, WorldPoint


def run(hydro, seed, dt):
    wec = WEC("OSWEC published held yaw, physical implicit path")
    flap = wec.body(
        "flap", hydro, mass=12700, inertia=(1.85e6,) * 3,
        passive_yaw=True, passive_yaw_threshold=1,
    )
    wec.body("base", hydro, mass=999, inertia=(999,) * 3)
    yaw = wec.coordinate("yaw", flap.move("yaw", pivot=WorldPoint(0, 0, -8.9)))
    wec.rotational_pto("hinge", yaw, damping=120000)
    result = wec.run(
        PMWave(2.5, 8, direction=10, seed=seed, phase_generator="matlab"),
        dt=dt, end_time=250, ramp_time=100, radiation_memory=40,
    )
    data = {
        "time": result.time,
        "wave": result.wave_elevation,
        "yaw": result.coordinates["yaw"].position,
        "speed": result.coordinates["yaw"].velocity,
        "torque": result.ptos["hinge"].force,
        "power": result.ptos["hinge"].absorbed_power,
    }
    expected_size = round(250 / dt) + 1
    assert all(len(value) == expected_size for value in data.values())
    assert all(np.isfinite(value).all() for value in data.values())
    np.testing.assert_allclose(data["torque"], -120000 * data["speed"],
                               rtol=0, atol=1e-7)
    assert np.trapezoid(data["power"], data["time"]) > 0
    return data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hydro", required=True, type=Path)
    parser.add_argument("--seed", required=True, type=int, choices=(1, 2, 3))
    parser.add_argument("--dt", required=True, type=float)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    data = run(args.hydro, args.seed, args.dt)
    np.savez_compressed(args.output, **data)
    print(f"Seed {args.seed}, dt={args.dt:g}: "
          f"physical PTO work {np.trapezoid(data['power'], data['time']):.6f} J")


if __name__ == "__main__":
    main()

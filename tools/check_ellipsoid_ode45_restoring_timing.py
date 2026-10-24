"""Check the causal force timing in a disposable ellipsoid source-model run.

This experiment changes only sampled subsystems in a copied MATLAB model.
It does not change or qualify the published source input.
"""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import trapezoid

from wecsim import RegularCICWave, RegularWave, WEC, WorldPoint
from wecsim.hydrodynamics.nonlinearHydro import HeaveMeshHydro


def _read(folder, mode, suffix, columns, rows):
    array = np.loadtxt(folder / f"{mode}_{suffix}.csv", delimiter=",", ndmin=2)
    assert array.shape == (rows, columns), (mode, suffix, array.shape)
    assert np.isfinite(array).all(), (mode, suffix)
    return array


def _max_error(a, b):
    assert a.shape == b.shape
    return float(np.max(np.abs(a - b)))


def _python_physical(wave_case, app):
    wec = WEC("ellipsoid force-timing source experiment")
    body = wec.body(
        "ellipsoid", app / "hydroData/ellipsoid.h5",
        mass="equilibrium", inertia=(1.375264e6, 1.375264e6, 1.341721e6),
        geometry_file=app / "geometry/elipsoid.stl",
        nonlinear_hydro="instantaneous", drag_coefficient=1,
        drag_area=np.pi * 25,
    )
    wec.coordinate("heave", body.move("heave"))
    wec.pto("PTO1", WorldPoint(0, 0, -12.5), body.at(0, 0, 0),
            damping=1_200_000)
    wave = (RegularCICWave(4, 6) if wave_case == "RegularCIC"
            else RegularWave(4, 6))
    return wec.run(
        wave, dt=.05, end_time=150, ramp_time=50, rho=1025,
        radiation_memory=60 if wave_case == "RegularCIC" else None,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("wave", choices=("Regular", "RegularCIC"))
    parser.add_argument("folder", type=Path)
    parser.add_argument("--duration", choices=("smoke", "full"), default="smoke")
    parser.add_argument("--applications", type=Path, default=Path("applications"))
    parser.add_argument("--modes", nargs="+", default=(
        "baseline", "continuous_restoring", "continuous_nonlinear"))
    args = parser.parse_args()
    stop = 2 if args.duration == "smoke" else 150
    rows = round(stop / .05) + 1
    geometry = args.applications / "Nonlinear_Hydro/geometry/elipsoid.stl"
    mesh = HeaveMeshHydro.from_stl(
        geometry, center_z=-2, rho=1025, gravity=9.81,
        depth=70, period=6, height=4, ramp_time=50,
        mass=None, drag_coefficient=1, drag_area=np.pi * 25,
    )
    result = {}
    for mode in args.modes:
        body = _read(args.folder, mode, "body", 55, rows)
        pto = _read(args.folder, mode, "pto", 25, rows)
        wave = _read(args.folder, mode, "wave", 2, rows)
        settings = _read(args.folder, mode, "settings", 9, 1)[0]
        expected_periods = {
            "baseline": (.05, .05, .05),
            "continuous_restoring": (0, .05, .05),
            "continuous_nonlinear": (0, 0, 0),
        }[mode]
        np.testing.assert_allclose(
            settings, [.05, .05, .05, .05, .05, stop, *expected_periods],
            rtol=0, atol=1e-12,
        )
        times = np.arange(rows) * .05
        for name, source in (("body", body), ("pto", pto), ("wave", wave)):
            assert _max_error(source[:, 0], times) < 1e-9, name
        heave = body[:, 3]
        speed = body[:, 9]
        source_restoring = -body[:, 39]
        current = np.array([
            mesh.forces(t, z + 2, v)[0]
            for t, z, v in zip(times, heave, speed, strict=True)
        ])
        previous = np.array([
            mesh.forces(t, z + 2, v)[0]
            for t, z, v in zip(times[:-1], heave[:-1], speed[:-1], strict=True)
        ])
        result[mode] = {
            "compiled_sample_periods_s": tuple(float(x) for x in settings[-3:]),
            "restoring_vs_current_N": _max_error(source_restoring[1:], current[1:]),
            "restoring_vs_previous_N": _max_error(source_restoring[1:], previous),
            "max_heave_m": float(np.max(np.abs(heave))),
            "max_speed_mps": float(np.max(np.abs(speed))),
            "wave": wave[:, 1],
            "heave": heave,
            "speed": speed,
            "pto_force": pto[:, 15],
            "pto_power": -pto[:, 15] * pto[:, 9],
        }
    baseline = result["baseline"]
    assert baseline["restoring_vs_previous_N"] < 1e-6, (
        "The copied baseline no longer has the archived preceding-state law"
    )
    assert baseline["restoring_vs_current_N"] > 5, (
        "The 2 s baseline is too short to distinguish the force clocks"
    )
    archive = (Path(__file__).resolve().parents[1]
               / "tests/fixtures/matlab/ellipsoid_ode45_nonlinear_clock"
               / "source_maxstep_common_grid.npz")
    prefix = "regular_cic" if args.wave == "RegularCIC" else "regular"
    with np.load(archive, allow_pickle=False) as original:
        for key, limit in (("heave", 1e-7), ("speed", 1e-7),
                           ("pto_force", .1), ("wave", 1e-10)):
            error = _max_error(
                baseline[key], original[f"{prefix}_published_{key}"][:rows])
            assert error < limit, ("baseline published trace", key, error)
    summary = {}
    if args.duration == "full":
        python = _python_physical(args.wave,
                                  args.applications / "Nonlinear_Hydro")
        python_body = python.bodies["ellipsoid"]
        python_pto = python.ptos["PTO1"]
        physical = {
            "heave": python_body.position[:, 2],
            "speed": python_body.velocity[:, 2],
            "pto_force": python_pto.force,
            "pto_power": python_pto.absorbed_power,
            "wave": python.wave_elevation,
        }
        assert _max_error(physical["wave"], baseline["wave"]) < 1e-10
    for mode, data in result.items():
        summary[mode] = {key: value for key, value in data.items()
                         if not isinstance(value, np.ndarray)}
        if mode != "baseline":
            summary[mode].update({
                "vs_baseline_heave_m": _max_error(data["heave"], baseline["heave"]),
                "vs_baseline_speed_mps": _max_error(data["speed"], baseline["speed"]),
                "vs_baseline_pto_force_N": _max_error(
                    data["pto_force"], baseline["pto_force"]),
                "vs_baseline_pto_power_W": _max_error(
                    data["pto_power"], baseline["pto_power"]),
                "vs_baseline_wave_m": _max_error(data["wave"], baseline["wave"]),
            })
            assert summary[mode]["vs_baseline_wave_m"] < 1e-10, (
                "The incident wave changed in a force-timing arm", mode
            )
        if args.duration == "full":
            summary[mode].update({
                f"vs_python_{key}": _max_error(data[key], physical[key])
                for key in ("heave", "speed", "pto_force", "pto_power")
            })
            source_work = float(trapezoid(data["pto_power"], times))
            python_work = float(trapezoid(physical["pto_power"], times))
            summary[mode]["source_signed_work_J"] = source_work
            summary[mode]["vs_python_signed_work_fraction"] = abs(
                source_work / python_work - 1)
            motion_limits = ((.0035, .0045, 6000, 5000)
                             if args.wave == "RegularCIC"
                             else (.003, .004, 5000, 4000))
            summary[mode]["passes_inherited_motion_and_work_limits"] = (
                all(summary[mode][f"vs_python_{key}"] < limit
                    for key, limit in zip(
                        ("heave", "speed", "pto_force", "pto_power"),
                        motion_limits, strict=True))
                and summary[mode]["vs_python_signed_work_fraction"] < .001
            )
    print(json.dumps({"wave": args.wave, "duration_s": stop,
                      "results": summary}, indent=2))


if __name__ == "__main__":
    main()

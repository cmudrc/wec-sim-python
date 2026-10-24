"""Screen copied-model continuous nonlinear ellipsoid source over ode45 steps.

The source input, 0.05 s output, nonlinear, and CIC clocks are held fixed.
Only the copied model's ode45 MaxStep changes from 0.05 to 0.025 to 0.0125 s.
This derived screen does not qualify the unchanged published input.
"""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import trapezoid

from wecsim import RegularCICWave, RegularWave, WEC, WorldPoint
from wecsim.hydrodynamics.nonlinearHydro import HeaveMeshHydro


SOURCE_LIMITS = {
    "Regular": (.001, .0013, 1600, 1300),
    "RegularCIC": (.0011, .0015, 2000, 1600),
}
PHYSICAL_LIMITS = {
    "Regular": (.003, .004, 5000, 4000),
    "RegularCIC": (.0035, .0045, 6000, 5000),
}
CHANNELS = ("heave", "speed", "pto_force", "pto_power")


def _read(folder, name, columns):
    values = np.loadtxt(folder / name, delimiter=",", ndmin=2)
    assert values.shape == (3001, columns), (folder, name, values.shape)
    assert np.isfinite(values).all(), (folder, name)
    return values


def _max_error(left, right):
    assert left.shape == right.shape
    return float(np.max(np.abs(left - right)))


def _source(folder, step, mesh):
    body = _read(folder, "continuous_nonlinear_body.csv", 55)
    pto = _read(folder, "continuous_nonlinear_pto.csv", 25)
    wave = _read(folder, "continuous_nonlinear_wave.csv", 2)
    settings = np.loadtxt(
        folder / "continuous_nonlinear_settings.csv", delimiter=",", ndmin=2)
    np.testing.assert_allclose(
        settings, [[.05, .05, .05, .05, step, 150, 0, 0, 0]],
        rtol=0, atol=1e-12,
    )
    time = np.arange(3001) * .05
    for array in (body, pto, wave):
        assert _max_error(array[:, 0], time) < 1e-9
    heave, speed, restoring = body[:, 3], body[:, 9], -body[:, 39]
    expected = np.array([
        mesh.forces(t, z + 2, v)[0]
        for t, z, v in zip(time, heave, speed, strict=True)
    ])
    assert _max_error(restoring, expected) < 1e-6, (
        "Copied-model restoring force no longer follows current state", folder
    )
    assert _max_error(pto[:, 9], speed) < 1e-9
    assert _max_error(pto[:, 15], -1_200_000 * pto[:, 9]) < 1e-6
    return {
        "time": time,
        "heave": heave,
        "speed": speed,
        "pto_force": pto[:, 15],
        "pto_power": -pto[:, 15] * pto[:, 9],
        "wave": wave[:, 1],
    }


def _python(case, app):
    wec = WEC("ellipsoid continuous nonlinear source-step screen")
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
    wave = RegularCICWave(4, 6) if case == "RegularCIC" else RegularWave(4, 6)
    output = wec.run(
        wave, dt=.05, end_time=150, ramp_time=50, rho=1025,
        radiation_memory=60 if case == "RegularCIC" else None,
    )
    return {
        "time": output.time,
        "heave": output.bodies["ellipsoid"].position[:, 2],
        "speed": output.bodies["ellipsoid"].velocity[:, 2],
        "pto_force": output.ptos["PTO1"].force,
        "pto_power": output.ptos["PTO1"].absorbed_power,
        "wave": output.wave_elevation,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", choices=("Regular", "RegularCIC"))
    parser.add_argument("base", type=Path)
    parser.add_argument("step025", type=Path)
    parser.add_argument("step0125", type=Path)
    parser.add_argument("--applications", type=Path, default=Path("applications"))
    args = parser.parse_args()
    app = args.applications / "Nonlinear_Hydro"
    mesh = HeaveMeshHydro.from_stl(
        app / "geometry/elipsoid.stl", center_z=-2, rho=1025,
        gravity=9.81, depth=70, period=6, height=4, ramp_time=50,
        mass=None, drag_coefficient=1, drag_area=np.pi * 25,
    )
    labels = ("maxstep_0p05", "maxstep_0p025", "maxstep_0p0125")
    source = {
        label: _source(folder, step, mesh)
        for label, folder, step in zip(
            labels, (args.base, args.step025, args.step0125),
            (.05, .025, .0125), strict=True
        )
    }
    reference = _python(args.case, app)
    first, middle, final = (source[label] for label in labels)
    for history in (first, middle, final):
        assert _max_error(history["wave"], first["wave"]) < 1e-10
        assert _max_error(history["wave"], reference["wave"]) < 1e-10
    summary = {"case": args.case, "source_step": {}, "physical": {}}
    source_ok = True
    physical_ok = True
    for name, source_limit, physical_limit in zip(
        CHANNELS, SOURCE_LIMITS[args.case], PHYSICAL_LIMITS[args.case],
        strict=True,
    ):
        coarse = _max_error(first[name], middle[name])
        fine = _max_error(middle[name], final[name])
        physical = _max_error(final[name], reference[name])
        summary["source_step"][name] = {
            "0p05_to_0p025": coarse,
            "0p025_to_0p0125": fine,
            "final_limit": source_limit,
            "contracts": fine < .7 * coarse,
            "passes_final_limit": fine < source_limit,
        }
        summary["physical"][name] = {
            "final_error": physical,
            "limit": physical_limit,
            "passes": physical < physical_limit,
        }
        source_ok &= fine < source_limit and fine < .7 * coarse
        physical_ok &= physical < physical_limit
    source_work = [float(trapezoid(history["pto_power"], history["time"]))
                   for history in (first, middle, final)]
    python_work = float(trapezoid(reference["pto_power"], reference["time"]))
    work_change = [abs(later - earlier) / later
                   for earlier, later in zip(source_work, source_work[1:])]
    work_error = abs(source_work[-1] / python_work - 1)
    summary["source_step"]["signed_work"] = {
        "work_J": source_work,
        "relative_changes": work_change,
        "passes_final": work_change[-1] < .001,
        "contracts": work_change[-1] < .7 * work_change[0],
    }
    summary["physical"]["signed_work"] = {
        "final_error_fraction": work_error,
        "limit": .001,
        "passes": work_error < .001,
    }
    source_ok &= work_change[-1] < .001 and work_change[-1] < .7 * work_change[0]
    physical_ok &= work_error < .001
    summary["passes_source_step_screen"] = bool(source_ok)
    summary["passes_inherited_physical_limits"] = bool(physical_ok)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

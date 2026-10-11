"""Full-duration refined ode45 nonlinear-hydro source and Python pair."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim import RegularCICWave, RegularWave, WEC, WorldPoint
from wecsim.nonlinearHydro import HeaveMeshHydro


APPLICATIONS = os.environ.get("WEC_SIM_APPLICATIONS_DIR")
REFERENCE = os.environ.get("WEC_SIM_ELLIPSOID_ODE45_FINE_DIR")
WAVE_CASE = os.environ.get("WEC_SIM_ELLIPSOID_CASE")
pytestmark = pytest.mark.skipif(
    not (APPLICATIONS and REFERENCE and WAVE_CASE),
    reason="pinned ellipsoid ode45 refinement output absent",
)

# Use half the established ode4 motion gates to check source refinement.
SOURCE_LIMITS = {
    "Regular": (.003, .00325, 4000, 4000),
    "RegularCIC": (.00375, .00425, 5000, 5000),
}
PYTHON_LIMITS = {
    "Regular": (.003, .004, 5000, 4000),
    "RegularCIC": (.0035, .0045, 6000, 5000),
}
PYTHON_STEP_LIMITS = (.0001, .0001, 100, 100)


def _read(step, kind):
    values = np.loadtxt(Path(REFERENCE) / WAVE_CASE / f"{step}_{kind}.csv",
                        delimiter=",", ndmin=2)
    assert np.isfinite(values).all(), f"{step} {kind}"
    return values


def _check(actual, expected, limit, label, failures):
    assert actual.shape == expected.shape, label
    error = float(np.max(np.abs(actual - expected)))
    print(f"{label}: {error:.8g} (gate {limit:.8g})")
    if not np.isfinite(error) or error >= limit:
        failures.append(f"{label}: {error:.8g} >= {limit:.8g}")


def _motion(actual_position, actual_speed, actual_force, actual_power,
            source_body, source_pto, limits, label, failures):
    _check(actual_position, source_body[:, 3], limits[0],
           f"{label} heave position", failures)
    _check(actual_speed, source_body[:, 9], limits[1],
           f"{label} heave speed", failures)
    _check(actual_force, source_pto[:, 15], limits[2],
           f"{label} PTO force", failures)
    _check(actual_power, -source_pto[:, 15] * source_pto[:, 9], limits[3],
           f"{label} absorbed power", failures)


def _configured_wec(app, geometry):
    wec = WEC("refined ode45 ellipsoid")
    ellipsoid = wec.body(
        "ellipsoid", app / "hydroData/ellipsoid.h5",
        mass="equilibrium", inertia=(1.375264e6, 1.375264e6, 1.341721e6),
        geometry_file=geometry, nonlinear_hydro="instantaneous",
        drag_coefficient=1, drag_area=np.pi * 25,
    )
    wec.coordinate("heave", ellipsoid.move("heave"))
    wec.pto("PTO1", WorldPoint(0, 0, -12.5), ellipsoid.at(0, 0, 0),
            damping=1_200_000)
    return wec


def test_refined_ode45_ellipsoid_full_trajectory():
    assert WAVE_CASE in SOURCE_LIMITS
    coarse_body, fine_body = (_read(step, "body")
                              for step in ("step010", "step005"))
    coarse_pto, fine_pto = (_read(step, "pto")
                            for step in ("step010", "step005"))
    coarse_wave, fine_wave = (_read(step, "wave")
                              for step in ("step010", "step005"))
    assert coarse_body.shape == (15001, 55)
    assert fine_body.shape == (30001, 55)
    assert coarse_pto.shape == (15001, 25)
    assert fine_pto.shape == (30001, 25)
    assert coarse_wave.shape == (15001, 2)
    assert fine_wave.shape == (30001, 2)
    failures = []
    _check(coarse_wave, fine_wave[::2], 1e-10,
           "source fine-step wave", failures)
    _motion(coarse_body[:, 3], coarse_body[:, 9], coarse_pto[:, 15],
            -coarse_pto[:, 15] * coarse_pto[:, 9],
            fine_body[::2], fine_pto[::2], SOURCE_LIMITS[WAVE_CASE],
            "source 0.01 to 0.005 s step change", failures)

    app = Path(APPLICATIONS) / "Nonlinear_Hydro"
    geometry = app / "geometry/elipsoid.stl"
    mesh = HeaveMeshHydro.from_stl(
        geometry, center_z=-2, rho=1025, gravity=9.81,
        depth=70, period=6, height=4, ramp_time=50,
        mass=None, drag_coefficient=1, drag_area=np.pi * 25,
    )
    for step, body in (("step010", coarse_body), ("step005", fine_body)):
        samples = np.arange(1, len(body), 10)
        previous = np.array([
            mesh.forces(body[i - 1, 0], body[i - 1, 3] + 2,
                        body[i - 1, 9])[0]
            for i in samples
        ])
        current = np.array([
            mesh.forces(body[i, 0], body[i, 3] + 2, body[i, 9])[0]
            for i in samples
        ])
        _check(previous, -body[samples, 39], 1e-6,
               f"{step} preceding-state restoring force", failures)
        print(f"{step} current-state restoring-force mismatch: "
              f"{np.max(np.abs(current + body[samples, 39])):.8g} N")

    def advance(dt):
        wave = (RegularCICWave(4, 6) if WAVE_CASE == "RegularCIC"
                else RegularWave(4, 6))
        return _configured_wec(app, geometry).run(
            wave, dt=dt, end_time=150, ramp_time=50,
            radiation_memory=60 if WAVE_CASE == "RegularCIC" else None,
            rho=1025,
        )

    coarse_python = advance(.01)
    fine_python = advance(.005)
    assert coarse_python.time.shape == (15001,)
    assert fine_python.time.shape == (30001,)
    _check(coarse_python.time, fine_python.time[::2], 1e-9,
           "Python common times", failures)
    for label, coarse_values, fine_values, limit in (
        ("heave position", coarse_python.bodies["ellipsoid"].position[:, 2],
         fine_python.bodies["ellipsoid"].position[::2, 2],
         PYTHON_STEP_LIMITS[0]),
        ("heave speed", coarse_python.bodies["ellipsoid"].velocity[:, 2],
         fine_python.bodies["ellipsoid"].velocity[::2, 2],
         PYTHON_STEP_LIMITS[1]),
        ("PTO force", coarse_python.ptos["PTO1"].force,
         fine_python.ptos["PTO1"].force[::2], PYTHON_STEP_LIMITS[2]),
        ("absorbed power", coarse_python.ptos["PTO1"].absorbed_power,
         fine_python.ptos["PTO1"].absorbed_power[::2],
         PYTHON_STEP_LIMITS[3]),
    ):
        _check(coarse_values, fine_values, limit,
               f"Python step change {label}", failures)
    _check(fine_python.time, fine_body[:, 0], 1e-9,
           "Python and source time", failures)
    _check(fine_python.wave_elevation, fine_wave[:, 1], 1e-10,
           "Python and source wave elevation", failures)
    _motion(fine_python.bodies["ellipsoid"].position[:, 2],
            fine_python.bodies["ellipsoid"].velocity[:, 2],
            fine_python.ptos["PTO1"].force,
            fine_python.ptos["PTO1"].absorbed_power,
            fine_body, fine_pto, PYTHON_LIMITS[WAVE_CASE],
            "Python versus refined source", failures)
    assert not failures, "; ".join(failures)

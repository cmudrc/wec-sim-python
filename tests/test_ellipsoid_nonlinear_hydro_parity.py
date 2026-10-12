"""Paired pinned Nonlinear_Hydro heaving-ellipsoid validation."""

import os
from pathlib import Path

import numpy as np
import pytest
from scipy.signal import fftconvolve

from wecsim import RegularCICWave, RegularWave, WEC, WorldPoint
from wecsim.hydrodynamics.bodyClass import BodyClass
from wecsim.hydrodynamics.nonlinearHydro import HeaveMeshHydro


APPLICATIONS = os.environ.get("WEC_SIM_APPLICATIONS_DIR")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
MODEL = os.environ.get("WEC_SIM_REFERENCE_MODEL", "ELLIPSOID_NLH_REG")
CASES = ({"ELLIPSOID_NLH_REG": ("ode4_Regular",),
          "ELLIPSOID_NLH_CIC": ("ode4_RegularCIC",),
          "ELLIPSOID_NLH_ODE45": ("ode45_Regular", "ode45_RegularCIC")}
         .get(MODEL, ()))
pytestmark = pytest.mark.skipif(
    not (APPLICATIONS and REFERENCE),
    reason="paired MATLAB nonlinear-hydro output and Applications absent",
)


def _source(case_name):
    source = Path(REFERENCE)
    body = np.loadtxt(source / f"{MODEL}_{case_name}_body1.csv",
                      delimiter=",")
    pto = np.loadtxt(source / f"{MODEL}_{case_name}_pto1.csv",
                     delimiter=",")
    mass_prefix = (f"{MODEL}_{case_name}" if MODEL == "ELLIPSOID_NLH_ODE45"
                   else MODEL)
    mass = np.loadtxt(source / f"{mass_prefix}_mass.csv", delimiter=",")
    assert body.shape == (3001, 55) and pto.shape == (3001, 25)
    return body, pto, mass


def _files():
    app = Path(APPLICATIONS) / "Nonlinear_Hydro"
    return app / "hydroData/ellipsoid.h5", app / "geometry/elipsoid.stl"


def _max_error(actual, expected, limit, label):
    assert actual.shape == expected.shape, label
    assert np.isfinite(actual).all() and np.isfinite(expected).all(), label
    error = np.max(np.abs(actual - expected))
    assert error < limit, f"{label}: {error:.6g} exceeds {limit}"


@pytest.mark.parametrize("case_name", CASES)
def test_mesh_forces_on_matlab_trajectory(case_name):
    body, pto, mass = _source(case_name)
    hydro_file, geometry_file = _files()
    mesh = HeaveMeshHydro.from_stl(
        geometry_file, center_z=-2, rho=1025, gravity=9.81,
        depth=70, period=6, height=4, ramp_time=50,
        mass=None, drag_coefficient=1, drag_area=np.pi * 25,
    )
    _max_error(np.array([mesh.mass]), mass[:1], 1e-6, "mesh equilibrium mass")
    forces = np.array([
        mesh.forces(t, z + 2, v)
        for t, z, v in zip(body[:, 0], body[:, 3], body[:, 9])
    ])
    if MODEL == "ELLIPSOID_NLH_ODE45":
        # The source's ode45 nonlinear buoyancy block holds the preceding
        # 0.05 s sample while motion and the other force logs advance.
        _max_error(forces[:-1, 0], -body[1:, 39], 1e-6,
                   "sampled mesh buoyancy minus weight")
        assert np.max(np.abs(forces[:, 0] + body[:, 39])) > 40_000
    else:
        _max_error(forces[:, 0], -body[:, 39], 1e-6,
                   "mesh buoyancy minus weight")
    _max_error(forces[:, 2], -body[:, 45], 1e-6, "quadratic drag")

    # WEC-Sim logs its linear BEM excitation plus the nonlinear FK correction.
    bem = BodyClass(str(hydro_file))
    bem.bodyNumber = bem.bodyTotal = 1
    bem.readH5file()
    bem.mass = mesh.mass
    bem.inertia = mass[1:].tolist()
    bem.hydroStiffness = np.zeros((6, 6))
    bem.viscDrag = {"Drag": np.zeros((6, 6)), "cd": np.zeros(6),
                    "characteristicArea": np.zeros(6)}
    bem.linearDamping = np.zeros((6, 6))
    omega = 2 * np.pi / 6
    cic = case_name.endswith("RegularCIC")
    convolution_time = np.arange(1201) * .05 if cic else np.array([0.0])
    bem.hydroForcePre(omega, [0], len(convolution_time), convolution_time,
                      [], .05, 1025, 9.81,
                      "regularCIC" if cic else "regular",
                      np.zeros((2, 2)), 1, 1, 0, 0, 0)
    re = np.asarray(bem.hydroForce["fExt"]["re"])[2]
    im = np.asarray(bem.hydroForce["fExt"]["im"])[2]
    ramp = np.array([mesh.ramp(t) for t in body[:, 0]])
    linear = 2 * ramp * (re * np.cos(omega * body[:, 0])
                         - im * np.sin(omega * body[:, 0]))
    _max_error(linear + forces[:, 1], body[:, 21], 1e-5,
               "linear plus nonlinear Froude-Krylov excitation")
    if MODEL == "ELLIPSOID_NLH_ODE45":
        # The lagged restoring force is in the applied source force sum,
        # not merely a reporting offset in the body output.
        total = (body[:, 21] - body[:, 27] - body[:, 33]
                 - body[:, 39] - body[:, 45])
        _max_error(total, body[:, 15], 1e-6, "source force assembly")
        adjusted_mass = (mesh.mass + 2 * np.trace(
            np.asarray(bem.hydroForce["fAddedMass"])[:3, :3]
        ))
        _max_error(adjusted_mass * body[:, 51], body[:, 15] + pto[:, 15],
                   1e-6, "source acceleration and adjusted mass")
    if cic:
        kernel = np.asarray(bem.hydroForce["irkb"])[:, 2, 2]
        velocity = body[:, 9]
        radiation = .05 * fftconvolve(velocity, kernel)[:len(velocity)]
        radiation -= .05 / 2 * kernel[0] * velocity
        last_lag = np.minimum(np.arange(len(velocity)), len(kernel) - 1)
        radiation -= .05 / 2 * kernel[last_lag] * (
            velocity[np.arange(len(velocity)) - last_lag]
        )
        _max_error(radiation, body[:, 27],
                   300 if MODEL == "ELLIPSOID_NLH_ODE45" else 1e-6,
                   "regularCIC radiation convolution")


@pytest.mark.parametrize("case_name", CASES)
def test_public_python_configuration_against_matlab_motion_and_pto(case_name):
    body, pto, mass = _source(case_name)
    hydro_file, geometry_file = _files()
    wec = WEC("ellipsoid")
    ellipsoid = wec.body(
        "ellipsoid", hydro_file, mass="equilibrium", inertia=mass[1:],
        geometry_file=geometry_file, nonlinear_hydro="instantaneous",
        drag_coefficient=1, drag_area=np.pi * 25,
    )
    wec.coordinate("heave", ellipsoid.move("heave"))
    wec.pto("PTO1", WorldPoint(0, 0, -12.5), ellipsoid.at(0, 0, 0),
            damping=1_200_000)
    cic = case_name.endswith("RegularCIC")
    result = wec.run(RegularCICWave(4, 6) if cic else RegularWave(4, 6),
                     dt=.05, end_time=150, ramp_time=50,
                     radiation_memory=60 if cic else None, rho=1025)
    np.testing.assert_allclose(result.time, body[:, 0], rtol=0, atol=1e-10)
    if MODEL == "ELLIPSOID_NLH_ODE45":
        # The published MATLAB ode4 and ode45 runs themselves differ by up
        # to 16 mm and 35 kN with identical physical settings. This is a
        # solver-envelope check, not a claim of matching ode45 internals.
        limits = ((.023, .034, 40_000, 33_000) if cic
                  else (.022, .032, 38_000, 30_000))
    else:
        limits = ((.0075, .0085, 10_000, 10_000) if cic
                  else (.006, .0065, 8_000, 8_000))
    _max_error(result.bodies["ellipsoid"].position[:, 2], body[:, 3],
               limits[0], "heave position")
    _max_error(result.bodies["ellipsoid"].velocity[:, 2], body[:, 9],
               limits[1], "heave velocity")
    _max_error(result.ptos["PTO1"].stroke, pto[:, 3], limits[0],
               "PTO stroke")
    _max_error(result.ptos["PTO1"].force, pto[:, 15], limits[2],
               "PTO force")
    _max_error(result.ptos["PTO1"].absorbed_power,
               -pto[:, 15] * pto[:, 9], limits[3], "PTO absorbed power")

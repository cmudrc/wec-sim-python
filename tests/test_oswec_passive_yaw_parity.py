"""Pair both published OSWEC passive-yaw cases with pinned MATLAB output."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim import RegularWave, WEC, WorldPoint
from wecsim.hydrodynamics.bodyClass import BodyClass
from wecsim.hydrodynamics.passiveYaw import PassiveYawExcitation


APPLICATIONS = os.environ.get("WEC_SIM_APPLICATIONS_DIR")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
pytestmark = pytest.mark.skipif(
    not (APPLICATIONS and REFERENCE),
    reason="paired MATLAB passive-yaw output not provided",
)


def _max_error(actual, expected, limit, label):
    actual = np.asarray(actual)
    expected = np.asarray(expected)
    assert actual.shape == expected.shape, label
    assert np.isfinite(expected).all(), label
    error = np.max(np.abs(actual - expected))
    assert error < limit, f"{label}: {error:.6g} exceeds {limit}"


@pytest.mark.parametrize("tracking", [False, True])
def test_published_passive_yaw_against_matlab(tracking):
    label = "PassiveYawON" if tracking else "PassiveYawOFF"
    source = Path(REFERENCE)
    hydro = (Path(APPLICATIONS)
             / "_Common_Input_Files/OSWEC/hydroData/oswec.h5").resolve()
    flap = np.loadtxt(source / f"OSWEC_PASSIVE_YAW_{label}_body1.csv",
                      delimiter=",")
    base = np.loadtxt(source / f"OSWEC_PASSIVE_YAW_{label}_body2.csv",
                      delimiter=",")
    pto = np.loadtxt(source / f"OSWEC_PASSIVE_YAW_{label}_pto1.csv",
                     delimiter=",")
    wave = np.loadtxt(source / f"OSWEC_PASSIVE_YAW_{label}_wave.csv",
                      delimiter=",")
    wec = WEC(f"OSWEC {label}")
    moving = wec.body("flap", hydro, mass=12700,
                      inertia=(1.85e6,) * 3, passive_yaw=tracking)
    wec.body("base", hydro, mass=999, inertia=(999,) * 3)
    yaw = wec.coordinate(
        "yaw", moving.move("yaw", pivot=WorldPoint(0, 0, -8.9)),
    )
    wec.rotational_pto("hinge", yaw, damping=120000)
    result = wec.run(RegularWave(2.5, 8, 10), dt=0.01,
                     end_time=600, ramp_time=100)

    assert flap.shape == base.shape == (60001, 25)
    assert pto.shape == (60001, 25)
    assert wave.shape == (60001, 2)
    np.testing.assert_allclose(result.time, flap[:, 0], rtol=0, atol=1e-10)
    np.testing.assert_allclose(result.time, pto[:, 0], rtol=0, atol=1e-10)
    _max_error(result.wave_elevation, wave[:, 1], 1e-12, "wave elevation")
    angle_limit = 0.002 if tracking else 1e-8
    speed_limit = 0.00025 if tracking else 1e-8
    _max_error(result.bodies["flap"].position[:, 5], flap[:, 6],
               angle_limit, f"{label} flap yaw")
    _max_error(result.bodies["flap"].velocity[:, 5], flap[:, 12],
               speed_limit, f"{label} flap yaw speed")
    _max_error(result.bodies["flap"].position[:, :5], flap[:, 1:6],
               1e-10, f"{label} stationary flap positions")
    _max_error(result.bodies["flap"].velocity[:, :5], flap[:, 7:12],
               1e-10, f"{label} stationary flap velocities")
    _max_error(result.bodies["base"].position, base[:, 1:7],
               1e-10, f"{label} fixed base position")
    _max_error(result.bodies["base"].velocity, base[:, 7:13],
               1e-10, f"{label} fixed base velocity")

    hinge = result.ptos["hinge"]
    _max_error(hinge.stroke, pto[:, 5], angle_limit, f"{label} PTO angle")
    _max_error(hinge.velocity, pto[:, 11], speed_limit,
               f"{label} PTO angular speed")
    _max_error(hinge.force, pto[:, 17],
               35 if tracking else 1e-4, f"{label} PTO torque")
    _max_error(-hinge.absorbed_power, pto[:, 23],
               0.7 if tracking else 1e-4, f"{label} PTO power")
    _max_error(pto[:, 17], -120000 * pto[:, 11], 1e-6,
               f"{label} source PTO law")

    if tracking:
        # MATLAB holds heading coefficients for 0.01 degree of yaw change;
        # Python interpolates continuously. The source force gate measures
        # that numerical difference instead of baking it into the model.
        body = BodyClass(str(hydro))
        body.bodyNumber = 1
        body.bodyTotal = 2
        body.readH5file()
        source_force = PassiveYawExcitation.from_hydro_data(
            body.hydroData, omega=2 * np.pi / 8,
            incident_direction=10, amplitude=1.25,
            ramp_time=100, rho=1000, g=9.81,
        )
        force = np.stack([
            source_force.force(t, angle)
            for t, angle in zip(result.time, flap[:, 6])
        ])
        for dof, limit in enumerate((100, 100, 1, 150, 70, 500)):
            _max_error(force[:, dof], flap[:, 19 + dof], limit,
                       f"{label} excitation DOF {dof}")

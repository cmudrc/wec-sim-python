"""Paired MATLAB validation for the reactive simple direct-drive Sphere PTO."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim import RegularWave, SimpleDirectDrive, WEC, WorldPoint
from wecsim.hydrodynamics.bodyClass import BodyClass


SPHERE_H5 = os.environ.get("WEC_SIM_SPHERE_H5")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
pytestmark = pytest.mark.skipif(
    not (SPHERE_H5 and REFERENCE),
    reason="paired MATLAB direct-drive output and Sphere HDF5 absent",
)


def _data():
    folder = Path(REFERENCE)
    body = np.loadtxt(folder / "SPHERE_REACTIVE_DDPTO_ReactiveWithPTO_body1.csv",
                      delimiter=",")
    pto = np.loadtxt(folder / "SPHERE_REACTIVE_DDPTO_ReactiveWithPTO_pto1.csv",
                     delimiter=",")
    controller = np.loadtxt(folder / "SPHERE_REACTIVE_DDPTO_controller.csv",
                            delimiter=",")
    drive = np.loadtxt(folder / "SPHERE_REACTIVE_DDPTO_drive.csv",
                       delimiter=",")
    wave = np.loadtxt(folder / "SPHERE_REACTIVE_DDPTO_wave.csv",
                      delimiter=",")
    mass = np.loadtxt(folder / "SPHERE_REACTIVE_DDPTO_mass.csv",
                      delimiter=",", ndmin=2)
    assert body.shape == (20_001, 49)
    assert pto.shape == (20_001, 25)
    assert controller.shape == (20_001, 13)
    assert drive.shape == (20_001, 11)
    assert wave.shape == (20_001, 2)
    assert mass.shape == (1, 5)
    return body, pto, controller, drive, wave, mass


def _max_error(actual, expected, limit, label):
    assert actual.shape == expected.shape, label
    assert np.isfinite(actual).all() and np.isfinite(expected).all(), label
    error = float(np.max(np.abs(actual - expected)))
    assert error < limit, f"{label}: {error:.6g} exceeds {limit}"


def test_matlab_controller_generator_and_force_balance():
    body, pto, controller, drive, wave, mass = _data()
    time, z, speed, acceleration = (body[:, 0], body[:, 3], body[:, 9],
                                   body[:, 45])
    shaft_speed, mechanical_power, shaft_torque = drive[:, 1:4].T
    inertia, friction, generator = drive[:, 4:7].T
    current, voltage, loss, electrical_power = drive[:, 7:11].T
    _max_error(controller[:, 3],
               -3.8e5 * speed + 1.52e5 * (z - z[0]),
               1e-5, "source reactive controller")
    _max_error(shaft_speed, 100 * speed, 1e-10, "shaft speed")
    # Simulink's Derivative block reports the shaft-speed change over the
    # preceding output interval. Its value is not the instantaneous body
    # acceleration logged at the same endpoint.
    _max_error(inertia[1:],
               -2 * 100 * np.diff(speed) / np.diff(time), 1e-8,
               "source sampled shaft inertia torque")
    _max_error(friction, -shaft_speed, 1e-10, "shaft friction torque")
    _max_error(shaft_torque, inertia + friction + generator,
               1e-8, "shaft torque")
    _max_error(current, generator / 7.186, 1e-9, "generator current")
    _max_error(voltage, 7.186 * shaft_speed, 1e-9, "generator voltage")
    _max_error(loss, .483 * current**2, 1e-8, "winding loss")
    _max_error(electrical_power, voltage * current + loss,
               1e-8, "source electrical power convention")
    _max_error(mechanical_power, shaft_speed * shaft_torque,
               1e-8, "source mechanical power")
    _max_error(body[:, 21] - body[:, 27] - body[:, 33] - body[:, 39],
               body[:, 15], 1e-5, "source hydrodynamic force assembly")
    _max_error(mass[0, 4] * acceleration,
               body[:, 15] + 100 * shaft_torque,
               1e-4, "source mass and PTO force balance")
    _max_error(pto[:, 3], z - z[0], 1e-10, "source PTO stroke")
    _max_error(pto[:, 9], speed, 1e-10, "source PTO speed")
    _max_error(pto[:, 15], np.zeros_like(time), 1e-7,
               "zero native PTO force")

    hydro = BodyClass(str(Path(SPHERE_H5)))
    hydro.bodyNumber = hydro.bodyTotal = 1
    hydro.readH5file()
    hydro.mass = float(mass[0, 0])
    hydro.hydroStiffness = np.zeros((6, 6))
    hydro.viscDrag = {"Drag": np.zeros((6, 6)), "cd": np.zeros(6),
                      "characteristicArea": np.zeros(6)}
    hydro.linearDamping = np.zeros((6, 6))
    hydro.hydroForcePre(2*np.pi/9.6664, [0], 1, np.array([0.0]), [],
                        .01, 1000, 9.81, "regular", np.zeros((2, 1)),
                        1, 1, 0, 0, 0)
    force = hydro.hydroForce
    ramp = np.ones_like(time)
    early = time < 50
    ramp[early] = (1 - np.cos(np.pi * time[early] / 50)) / 2
    excitation = 1.25 * ramp * (
        force["fExt"]["re"][2] * np.cos(2*np.pi*time/9.6664)
        - force["fExt"]["im"][2] * np.sin(2*np.pi*time/9.6664)
    )
    _max_error(excitation, body[:, 21], 1e-5, "BEM heave excitation")
    _max_error(force["fDamping"][2, 2] * speed, body[:, 27],
               1e-5, "BEM radiation damping")
    _max_error(force["linearHydroRestCoef"][2, 2] * (z - z[0]),
               body[:, 39], 1e-5, "BEM restoring")
    _max_error(wave[:, 1], 1.25 * ramp * np.cos(2*np.pi*time/9.6664),
               1e-12, "regular wave")


def test_public_direct_drive_pto_tracks_matlab_dynamics_and_power():
    body, pto, controller, drive, wave, _ = _data()
    wec = WEC("reactive sphere with direct drive")
    sphere = wec.body(
        "sphere", Path(SPHERE_H5), mass="equilibrium",
        inertia=(20_907_301, 21_306_090.66, 37_085_481.11),
    )
    wec.coordinate("heave", sphere.move("heave"))
    wec.pto(
        "PTO1", WorldPoint(0, 0, -2), sphere.at(0, 0, 0),
        axis=(0, 0, 1),
        direct_drive=SimpleDirectDrive(
            kp=3.8e5, ki=-1.52e5, torque_constant=7.186,
            gear_ratio=100, drivetrain_inertia=2, drivetrain_friction=1,
            winding_resistance=.483, winding_inductance=5.223e-3,
        ),
    )
    result = wec.run(RegularWave(height=2.5, period=9.6664),
                     dt=.01, end_time=200, ramp_time=50)
    motion = result.bodies["sphere"]
    pto_result = result.ptos["PTO1"]
    dd = pto_result.direct_drive
    assert dd is not None
    _max_error(result.time, body[:, 0], 1e-10, "time")
    _max_error(motion.position[:, 2], body[:, 3], .0002,
               "sphere heave")
    _max_error(motion.velocity[:, 2], body[:, 9], .0001,
               "sphere heave velocity")
    _max_error(pto_result.stroke, pto[:, 3], .0002, "PTO stroke")
    _max_error(dict(result.raw.extra_outputs)["body1_controller_force"],
               controller[:, 3], 50, "reactive controller force")
    _max_error(dd.shaft_velocity, drive[:, 1], .01, "shaft speed")
    _max_error(dd.mechanical_power, drive[:, 2], 50,
               "mechanical shaft power")
    _max_error(dd.shaft_torque, drive[:, 3], 1, "shaft torque")
    _max_error(dd.generator_torque, drive[:, 6], 1,
               "generator torque")
    _max_error(dd.current, drive[:, 7], .1, "generator current")
    _max_error(dd.voltage, drive[:, 8], .1, "generator voltage")
    _max_error(dd.resistance_loss, drive[:, 9], 30, "winding loss")
    _max_error(dd.electrical_power, drive[:, 10], 50,
               "source-signed electrical power")
    _max_error(pto_result.force, 100 * drive[:, 3], 100,
               "PTO force on body")
    _max_error(pto_result.absorbed_power, -drive[:, 2], 50,
               "positive absorbed mechanical power")
    _max_error(result.wave_elevation, wave[:, 1], 1e-12,
               "wave elevation")

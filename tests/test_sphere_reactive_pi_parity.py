"""Compare the published Sphere reactive controller with paired MATLAB output.

The pinned Simulink controller applies -Kp * heave_velocity - Ki * heave_position
to the translational PTO, with Ki negative in this published example. The
public PTO builder keeps that active position gain separate from passive
stiffness.
"""

import os
from pathlib import Path

import numpy as np

from wecsim import RegularWave, WEC, WorldPoint


FIXTURES = Path(__file__).parent / "fixtures/matlab"
SPHERE_H5 = (os.environ.get("WEC_SIM_SPHERE_H5")
             or FIXTURES / "sphere_pto_config/sphere.h5")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")

KP = 49_181.0
KI = -573_350.0


def _max_error(actual, expected, limit, label):
    assert actual.shape == expected.shape, label
    assert np.isfinite(expected).all(), label
    error = np.max(np.abs(actual - expected))
    assert error < limit, f"{label}: max error {error:.6g} exceeds {limit}"


def _source():
    if REFERENCE:
        reference = Path(REFERENCE)
        body = np.loadtxt(
            reference / "Sphere_Reactive_PI_Reactive_PI_body1.csv", delimiter=",",
        )
        controller = np.loadtxt(
            reference / "Sphere_Reactive_PI_controller.csv", delimiter=",",
        )
        assert body.shape == (10_001, 25)
        assert controller.shape == (10_001, 13)
        np.testing.assert_allclose(body[:, 0], controller[:, 0], rtol=0, atol=1e-10)
        return (body[:, 0], body[:, 3], body[:, 9],
                controller[:, 3], controller[:, 9])
    with np.load(FIXTURES / "sphere_reactive_pi.npz") as fixture:
        return tuple(fixture[key] for key in (
            "time", "position_heave", "velocity_heave",
            "controller_force", "controller_power",
        ))


def test_published_sphere_reactive_pi_against_matlab():
    wec = WEC("Sphere reactive PI")
    sphere = wec.body(
        "sphere", Path(SPHERE_H5).resolve(),
        inertia=(20_907_301, 21_306_090.66, 37_085_481.11),
    )
    wec.coordinate("heave", sphere.move("heave"))
    wec.pto(
        "reactive", WorldPoint(0, 0, -2), sphere.at(0, 0, 0),
        axis=(0, 0, 1), velocity_control_gain=KP,
        position_control_gain=KI,
    )
    response = wec.run(
        RegularWave(2.5, 9.6664), dt=0.02, end_time=200, ramp_time=50,
    )
    time, position, velocity, source_force, source_power = _source()
    assert time.shape == (10_001,)
    np.testing.assert_allclose(response.time, time, rtol=0, atol=1e-10)

    displacement = response.coordinates["heave"].position
    speed = response.coordinates["heave"].velocity
    force = -KP * speed - KI * displacement
    power = force * speed
    pto = response.ptos["reactive"]
    np.testing.assert_allclose(pto.force, force, rtol=0, atol=1e-8)
    np.testing.assert_allclose(pto.control_force, force, rtol=0, atol=1e-8)
    np.testing.assert_allclose(pto.absorbed_power, -power, rtol=0, atol=1e-7)
    np.testing.assert_allclose(pto.control_absorbed_power, -power, rtol=0, atol=1e-7)
    np.testing.assert_allclose(
        response.raw.pto_generalized_force[:, 0], pto.force, rtol=0, atol=1e-8,
    )
    _max_error(response.bodies["sphere"].position[:, 2], position, 0.010,
               "Sphere heave position")
    _max_error(response.bodies["sphere"].velocity[:, 2], velocity, 0.006,
               "Sphere heave velocity")
    _max_error(force, source_force, 5_500, "reactive controller force")
    _max_error(power, source_power, 40_000, "reactive controller power")

    # The saved MATLAB controller must also follow its documented gains.
    matlab_force = -KP * velocity - KI * (position - position[0])
    _max_error(matlab_force, source_force, 1e-4,
               "MATLAB controller gain relation")
    _max_error(source_force * velocity, source_power, 1e-3,
               "MATLAB controller power relation")

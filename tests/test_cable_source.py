"""Paired force law for the published MBARI cable application."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim import PlanarCableAttachment, WecSimCableTension
from wecsim.hydrodynamics.bodyClass import BodyClass


REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
HYDRO = os.environ.get("WEC_SIM_CABLE_H5")
pytestmark = pytest.mark.skipif(
    not REFERENCE, reason="pinned MBARI cable MATLAB output absent",
)


def _load(name):
    return np.loadtxt(Path(REFERENCE) / f"CABLE_SOURCE_{name}.csv",
                      delimiter=",", ndmin=2)


def test_published_mbari_cable_force_from_relative_motion():
    parameters = _load("parameters").ravel()
    np.testing.assert_allclose(parameters, (1e6, 100, 17.8, 18),
                               rtol=0, atol=1e-10)
    cable = _load("cable1")
    wave = _load("wave")
    bodies = [_load(f"Cable_body{number}") for number in (1, 2, 3)]
    assert cable.shape == (12_001, 37)
    assert wave.shape == (12_001, 2)
    for body in bodies:
        assert body.shape == (12_001, 25)
        np.testing.assert_allclose(body[:, 0], cable[:, 0],
                                   rtol=0, atol=1e-10)

    time = cable[:, 0]
    ramp = np.where(time >= 40, 1,
                    (1 + np.cos(np.pi + np.pi * time / 40)) / 2)
    expected_wave = .5 * ramp * np.cos(2 * np.pi * time / 8)
    np.testing.assert_allclose(wave[:, 0], time, rtol=0, atol=1e-10)
    np.testing.assert_allclose(wave[:, 1], expected_wave,
                               rtol=0, atol=1e-11)

    law = WecSimCableTension(*parameters)
    force = law.force_z(cable[:, 3], cable[:, 9])
    source_force = cable[:, 21]
    assert np.max(np.abs(source_force)) > 1e5
    np.testing.assert_allclose(force, source_force, rtol=0, atol=1e-5)
    assert np.all(force[np.abs(cable[:, 3] + parameters[3])
                        <= parameters[2]] == 0)


def test_published_cable_attachment_from_body_local_points():
    cable = _load("cable1")
    base = _load("Cable_body2")
    follower = _load("Cable_body3")
    attachment = PlanarCableAttachment(
        base_offset=(0, 1.95), follower_offset=(0, -5.2),
        initial_length=18,
    )

    def pose_and_rate(body):
        return (body[:, (1, 3, 5)], body[:, (7, 9, 11)])

    displacement, speed = attachment.motion(
        *pose_and_rate(base), *pose_and_rate(follower),
    )
    np.testing.assert_allclose(displacement, cable[:, 3],
                               rtol=0, atol=1e-10)
    np.testing.assert_allclose(speed, cable[:, 9],
                               rtol=0, atol=1e-11)
    force = WecSimCableTension(1e6, 100, 17.8, 18).force_z(
        displacement, speed,
    )
    np.testing.assert_allclose(force, cable[:, 21], rtol=0, atol=1e-5)

    base_wrench, follower_wrench = attachment.world_wrenches(
        base[:, (1, 3, 5)], follower[:, (1, 3, 5)], force,
    )
    np.testing.assert_allclose(base_wrench[0], (0, 200_000, 0),
                               rtol=0, atol=1e-8)
    np.testing.assert_allclose(follower_wrench[0], (0, -200_000, 0),
                               rtol=0, atol=1e-8)
    # Axial power must equal work done on both translating and pitching bodies.
    power = (np.sum(base_wrench * base[:, (7, 9, 11)], axis=1)
             + np.sum(follower_wrench * follower[:, (7, 9, 11)], axis=1))
    np.testing.assert_allclose(power, force * speed, rtol=0, atol=1e-7)


def test_published_cable_endpoint_drag_and_axial_joint_loads():
    """The source cable has two 1 kg drag bodies, including when slack."""
    cable = _load("cable1")
    base, follower = (_load("Cable_body2"), _load("Cable_body3"))
    base_accel, follower_accel = (
        _load("body2_forces")[:, 31:37],
        _load("body3_forces")[:, 31:37],
    )
    joint2, joint3 = (_load("constraint2"), _load("constraint3"))
    attachment = PlanarCableAttachment(
        base_offset=(0, 1.95), follower_offset=(0, -5.2),
        initial_length=18,
    )
    base_wrench, follower_wrench = attachment.world_wrenches(
        base[:, (1, 3, 5)], follower[:, (1, 3, 5)], cable[:, 21],
    )

    def endpoint_vertical(body, acceleration, offset):
        angle, speed = body[:, 5], body[:, 11]
        velocity = body[:, 9] - offset * np.sin(angle) * speed
        end_acceleration = (acceleration[:, 2] - offset
                            * (np.sin(angle) * acceleration[:, 4]
                               + np.cos(angle) * speed**2))
        drag = -0.5 * 1000 * 1.4 * 10 * np.abs(velocity) * velocity
        return drag, end_acceleration

    base_drag, base_acceleration = endpoint_vertical(
        base, base_accel, 1.95,
    )
    follower_drag, follower_acceleration = endpoint_vertical(
        follower, follower_accel, -5.2,
    )
    # Source joint 2 reports the reaction on the base drag body; joint 3
    # reports the opposite convention on the follower drag body.
    lower_reaction = -base_wrench[:, 1] - base_drag + base_acceleration
    upper_reaction = (follower_wrench[:, 1] + follower_drag
                      - follower_acceleration)
    np.testing.assert_allclose(lower_reaction, joint2[:, 21],
                               rtol=0, atol=2)
    np.testing.assert_allclose(upper_reaction, joint3[:, 21],
                               rtol=0, atol=7)
    slack = cable[:, 21] == 0
    assert np.mean(slack) > 0.3
    assert np.max(np.abs(joint2[slack, 21])) > 1_000


@pytest.mark.skipif(not HYDRO, reason="pinned MBARI hydrodynamics absent")
def test_published_mbari_body_force_components():
    """Pair the source's force paths on prescribed three-body motion."""
    time = _load("Cable_body1")[:, 0]
    coefficients = (
        ((1.15, 1.15, 1, .5, .5, 0),
         (2.9568, 2.9568, 5.4739, 5.4739, 5.4739, 0)),
        ((.8, .8, 0, 1, 1, 0), (1.12, 1.12, 6.56, 6.56, 6.56, 0)),
        ((1.15, 1.15, 0, 1.15, 1.15, 0),
         (1.496, 1.496, 0, 1.496, 1.496, 0)),
    )
    for index, mass in enumerate((1080, 815, 600), start=1):
        source_body = _load(f"Cable_body{index}")
        source_force = _load(f"body{index}_forces")
        assert source_force.shape == (12_001, 37)
        np.testing.assert_allclose(source_force[:, 0], time,
                                   rtol=0, atol=1e-10)
        velocity = source_body[:, 7:13]
        cd, area = (np.asarray(value) for value in coefficients[index - 1])
        viscous = .5 * 1000 * cd * area * np.abs(velocity) * velocity
        np.testing.assert_allclose(viscous, source_force[:, 19:25],
                                   rtol=0, atol=1e-7)
        np.testing.assert_allclose(
            source_body[:, 19:25] - sum(source_force[:, start:start + 6]
                                         for start in (1, 7, 13, 19, 25)),
            source_body[:, 13:19], rtol=0, atol=1e-7,
        )
        if index == 3:
            restoring = np.zeros_like(velocity)
            restoring[:, 2] = (mass - 1000 * .2) * 9.81
            np.testing.assert_allclose(restoring, source_force[:, 13:19],
                                       rtol=0, atol=1e-8)
            continue

        body = BodyClass(HYDRO)
        body.bodyNumber = index
        body.bodyTotal = 2
        body.readH5file()
        body.mass = mass
        body.hydroStiffness = np.zeros((6, 6))
        body.viscDrag = {"Drag": np.zeros((6, 6)), "cd": np.zeros(6),
                         "characteristicArea": np.zeros(6)}
        body.linearDamping = np.zeros((6, 6))
        body.hydroForcePre(
            2 * np.pi / 8, [0], 1, np.array([0.]), [], .01, 1000, 9.81,
            "regular", np.vstack((time, np.zeros_like(time))),
            index, 2, 0, 0, 0,
        )
        hydro = body.hydroForce
        np.testing.assert_allclose(hydro["fAddedMass"],
                                   _load(f"body{index}_added_mass"),
                                   rtol=0, atol=1e-8)
        np.testing.assert_allclose(hydro["fDamping"],
                                   _load(f"body{index}_radiation_damping"),
                                   rtol=0, atol=1e-8)
        radiation = velocity @ hydro["fDamping"].T
        restoring = ((source_body[:, 1:7] - source_body[0, 1:7])
                     @ hydro["linearHydroRestCoef"].T)
        restoring[:, 2] -= (1000 * float(np.asarray(body.dispVol).item())
                             - mass) * 9.81
        np.testing.assert_allclose(radiation, source_force[:, 1:7],
                                   rtol=0, atol=1e-7)
        np.testing.assert_allclose(restoring, source_force[:, 13:19],
                                   rtol=0, atol=1e-7)

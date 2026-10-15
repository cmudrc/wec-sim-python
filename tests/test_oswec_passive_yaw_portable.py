"""Keep both published 600 s OSWEC regular passive-yaw gates in CI."""

from pathlib import Path

import numpy as np
import pytest

from wecsim import RegularWave, WEC, WorldPoint
from wecsim.hydrodynamics.bodyClass import BodyClass
from wecsim.hydrodynamics.passiveYaw import PassiveYawExcitation


ROOT = Path(__file__).resolve().parents[1]
HYDRO = ROOT / "tests/fixtures/matlab/oswec_multi_wave/oswec.h5"
SOURCE = Path(__file__).parent / "fixtures/matlab/oswec_passive_yaw"


def _max_error(actual, expected, limit, label):
    actual = np.asarray(actual)
    expected = np.asarray(expected)
    assert actual.shape == expected.shape, label
    assert np.isfinite(actual).all() and np.isfinite(expected).all(), label
    error = float(np.max(np.abs(actual - expected)))
    assert error < limit, f"{label}: {error:.6g} exceeds {limit}"


@pytest.mark.parametrize("tracking", [False, True])
def test_published_passive_yaw_from_pinned_matlab(tracking):
    label = "on" if tracking else "off"
    with np.load(SOURCE / f"{label}.npz") as saved:
        source = {name: saved[name] for name in saved.files}
    assert source["time"].shape == (60001,)
    assert all(np.isfinite(value).all() for value in source.values())

    wec = WEC(f"OSWEC PassiveYaw{label.upper()}")
    moving = wec.body("flap", HYDRO, mass=12700,
                      inertia=(1.85e6,) * 3, passive_yaw=tracking)
    wec.body("base", HYDRO, mass=999, inertia=(999,) * 3)
    yaw = wec.coordinate(
        "yaw", moving.move("yaw", pivot=WorldPoint(0, 0, -8.9)),
    )
    wec.rotational_pto("hinge", yaw, damping=120000)
    result = wec.run(RegularWave(2.5, 8, 10), dt=.01,
                     end_time=600, ramp_time=100)

    _max_error(result.time, source["time"], 1e-10, "time")
    _max_error(result.wave_elevation, source["wave_elevation"],
               1e-12, "wave elevation")
    angle_limit = .002 if tracking else 1e-8
    speed_limit = .00025 if tracking else 1e-8
    flap = result.bodies["flap"]
    _max_error(flap.position[:, 5], source["yaw"], angle_limit, "flap yaw")
    _max_error(flap.velocity[:, 5], source["yaw_speed"], speed_limit,
               "flap yaw speed")

    # The source's inactive axes are stationary to below 1e-15; retain
    # their initial pose and measured maximum drift instead of full histories.
    for name, actual in (
        ("flap_position", flap.position[:, :5]),
        ("flap_velocity", flap.velocity[:, :5]),
        ("base_position", result.bodies["base"].position),
        ("base_velocity", result.bodies["base"].velocity),
    ):
        assert source[f"{name}_source_max_deviation"] < 1e-10
        reference = np.broadcast_to(source[f"{name}_reference"], actual.shape)
        _max_error(actual, reference, 1e-10, name)

    hinge = result.ptos["hinge"]
    _max_error(hinge.stroke, source["pto_stroke"], angle_limit, "PTO angle")
    _max_error(hinge.velocity, source["pto_velocity"], speed_limit,
               "PTO angular speed")
    _max_error(hinge.force, source["pto_force"],
               35 if tracking else 1e-4, "PTO torque")
    _max_error(-hinge.absorbed_power, source["pto_power"],
               .7 if tracking else 1e-4, "PTO power")
    _max_error(source["pto_force"], -120000 * source["pto_velocity"],
               1e-6, "source PTO law")

    if tracking:
        body = BodyClass(str(HYDRO))
        body.bodyNumber = 1
        body.bodyTotal = 2
        body.readH5file()
        excitation = PassiveYawExcitation.from_hydro_data(
            body.hydroData, omega=2 * np.pi / 8,
            incident_direction=10, amplitude=1.25,
            ramp_time=100, rho=1000, g=9.81,
        )
        force = np.stack([
            excitation.force(t, angle)
            for t, angle in zip(result.time, source["yaw"])
        ])
        for dof, limit in enumerate((100, 100, 1, 150, 70, 500)):
            _max_error(force[:, dof], source["excitation"][:, dof],
                       limit, f"excitation DOF {dof}")

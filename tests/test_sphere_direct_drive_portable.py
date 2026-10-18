"""Keep the published Sphere reactive direct-drive trajectory in ordinary CI."""

from pathlib import Path

import numpy as np

from wecsim import RegularWave, SimpleDirectDrive, WEC, WorldPoint


FIXTURE = Path(__file__).parent / "fixtures/matlab/sphere_reactive_drive"


def _close(name, actual, source, limit):
    assert actual.shape == source.shape == (20_001,), name
    assert np.isfinite(actual).all() and np.isfinite(source).all(), name
    error = float(np.max(np.abs(actual - source)))
    assert error < limit, f"{name}: {error:.6g} exceeds {limit}"


def test_published_reactive_direct_drive_through_public_wec():
    with np.load(FIXTURE / "source.npz") as source:
        time = source["time"]
        assert time.shape == (20_001,)
        np.testing.assert_allclose(time, np.arange(20_001) * 0.01,
                                   rtol=0, atol=1e-8)
        assert np.max(np.abs(source["shaft_torque"])) > 0

        wec = WEC("Sphere reactive direct drive")
        sphere = wec.body(
            "sphere", FIXTURE / "sphere.h5", mass="equilibrium",
            inertia=(20_907_301, 21_306_090.66, 37_085_481.11),
        )
        wec.coordinate("heave", sphere.move("heave"))
        wec.pto(
            "PTO1", WorldPoint(0, 0, -2), sphere.at(0, 0, 0),
            axis=(0, 0, 1),
            direct_drive=SimpleDirectDrive(
                kp=3.8e5, ki=-1.52e5, torque_constant=7.186,
                gear_ratio=100, drivetrain_inertia=2,
                drivetrain_friction=1, winding_resistance=.483,
                winding_inductance=5.223e-3,
            ),
        )
        result = wec.run(RegularWave(2.5, 9.6664), dt=.01,
                         end_time=200, ramp_time=50)
        np.testing.assert_allclose(result.time, time, rtol=0, atol=1e-10)
        motion = result.bodies["sphere"]
        pto = result.ptos["PTO1"]
        drive = pto.direct_drive
        assert drive is not None

        _close("heave position", motion.position[:, 2],
               source["body_heave"], .0002)
        _close("heave speed", motion.velocity[:, 2],
               source["body_speed"], .0001)
        _close("PTO stroke", pto.stroke, source["pto_stroke"], .0002)
        _close("controller force",
               dict(result.raw.extra_outputs)["body1_controller_force"],
               source["controller_force"], 50)
        for field, limit in (
            ("shaft_velocity", .01), ("mechanical_power", 50),
            ("shaft_torque", 1), ("generator_torque", 1),
            ("current", .1), ("voltage", .1),
            ("resistance_loss", 30), ("electrical_power", 50),
        ):
            _close(field, getattr(drive, field), source[field], limit)
        _close("PTO force", pto.force, 100 * source["shaft_torque"], 100)
        _close("absorbed power", pto.absorbed_power,
               -source["mechanical_power"], 50)
        _close("wave elevation", result.wave_elevation,
               source["wave_elevation"], 1e-12)

        # The archived generator channels must obey their electrical units.
        np.testing.assert_allclose(
            source["current"] * 7.186, source["generator_torque"],
            rtol=0, atol=1e-9,
        )
        np.testing.assert_allclose(
            source["electrical_power"],
            source["voltage"] * source["current"]
            + source["resistance_loss"],
            rtol=0, atol=1e-8,
        )


def test_direct_drive_public_initial_state_reaches_pto_and_generator():
    """A configured heave release must reach the coupled drive, not just the body."""
    wec = WEC("Sphere direct-drive release")
    # Use the fixture's displaced-volume mass through the explicit mass route.
    sphere = wec.body("sphere", FIXTURE / "sphere.h5", mass=261_724)
    wec.coordinate("heave", sphere.move("heave"))
    wec.pto(
        "PTO1", WorldPoint(0, 0, -2), sphere.at(0, 0, 0),
        axis=(0, 0, 1),
        direct_drive=SimpleDirectDrive(
            kp=3.8e5, ki=-1.52e5, torque_constant=7.186,
            gear_ratio=100, drivetrain_inertia=2,
            drivetrain_friction=1, winding_resistance=.483,
            winding_inductance=5.223e-3,
        ),
    )
    result = wec.run(
        RegularWave(0, 9.6664), dt=.01, end_time=.2, ramp_time=0,
        initial_coordinate={"heave": .3}, initial_speed={"heave": -.2},
    )
    body = result.bodies["sphere"]
    pto = result.ptos["PTO1"]
    drive = pto.direct_drive
    assert drive is not None
    np.testing.assert_allclose(result.coordinates["heave"].position[0],
                               .3, atol=1e-12)
    np.testing.assert_allclose(body.velocity[0, 2], -.2, atol=1e-12)
    np.testing.assert_allclose(pto.stroke[0], .3, atol=1e-12)
    np.testing.assert_allclose(pto.velocity[0], -.2, atol=1e-12)
    np.testing.assert_allclose(drive.shaft_velocity[0], -20, atol=1e-12)
    np.testing.assert_allclose(drive.current[0], 0, atol=1e-12)
    np.testing.assert_allclose(drive.voltage[0], -143.72, atol=1e-12)
    np.testing.assert_allclose(
        dict(result.raw.extra_outputs)["body1_controller_force"][0],
        121_600, atol=1e-9,
    )
    np.testing.assert_allclose(
        pto.absorbed_power, -drive.mechanical_power, rtol=0, atol=1e-8,
    )

"""Short public-API checks for heave modes whose long MATLAB pairs run separately."""

from pathlib import Path

import numpy as np

from wecsim import HydroState, RegularWave, SimpleDirectDrive, WEC, WorldPoint
from wecsim.hydrodynamics.bodyClass import BodyClass


SPHERE = Path(__file__).resolve().parents[1] / "examples/data/sphere.h5"


def test_short_direct_drive_preserves_generator_and_pto_power_balance():
    wec = WEC("short direct-drive sphere")
    sphere = wec.body("sphere", SPHERE)
    wec.coordinate("heave", sphere.move("heave"))
    wec.pto(
        "drive", WorldPoint(0, 0, -2), sphere.at(0, 0, 0),
        axis=(0, 0, 1),
        direct_drive=SimpleDirectDrive(
            kp=3.8e5, ki=-1.52e5, torque_constant=7.186,
            gear_ratio=100, drivetrain_inertia=2, drivetrain_friction=1,
            winding_resistance=.483, winding_inductance=5.223e-3,
        ),
    )
    result = wec.run(RegularWave(2.5, 8), dt=.01, end_time=.2, ramp_time=0)
    drive = result.ptos["drive"].direct_drive
    assert drive is not None
    assert len(result.time) == 21
    assert np.isfinite(result.bodies["sphere"].position).all()
    np.testing.assert_allclose(
        drive.shaft_velocity, 100 * result.bodies["sphere"].velocity[:, 2],
    )
    np.testing.assert_allclose(drive.generator_torque, 7.186 * drive.current)
    np.testing.assert_allclose(drive.resistance_loss, .483 * drive.current**2)
    np.testing.assert_allclose(
        drive.mechanical_power, drive.shaft_velocity * drive.shaft_torque,
    )
    np.testing.assert_allclose(
        result.ptos["drive"].absorbed_power, -drive.mechanical_power,
    )


def _variable_sphere(second_mass, *, height=0, ramp_time=0,
                     switch_time=.1, end_time=.2):
    wec = WEC("two-state sphere")
    sphere = wec.variable_body(
        "sphere", (HydroState(SPHERE), HydroState(SPHERE, second_mass)),
        switch_times=(switch_time,),
    )
    wec.coordinate("heave", sphere.move("heave"))
    wec.pto("damper", WorldPoint(0, 0, 2), sphere.at(0, 0, 0),
            axis=(0, 0, 1), damping=200_000)
    return wec.run(RegularWave(height, 8), dt=.01,
                   end_time=end_time, ramp_time=ramp_time)


def test_short_variable_hydro_switch_changes_motion_only_after_switch():
    hydro = BodyClass(str(SPHERE))
    hydro.bodyNumber = hydro.bodyTotal = 1
    hydro.readH5file()
    equilibrium_mass = 1000 * float(np.asarray(hydro.dispVol).item())

    changed = _variable_sphere(1.1 * equilibrium_mass)
    control = _variable_sphere("equilibrium")
    active = dict(changed.raw.extra_outputs)["body1_hydro_state"]
    np.testing.assert_array_equal(active[:10], np.ones(10))
    np.testing.assert_array_equal(active[10:], np.full(11, 2))
    np.testing.assert_allclose(
        changed.bodies["sphere"].position[:10],
        control.bodies["sphere"].position[:10], rtol=0, atol=1e-12,
    )
    assert abs(changed.bodies["sphere"].position[-1, 2]
               - control.bodies["sphere"].position[-1, 2]) > 1e-8
    assert np.isfinite(changed.ptos["damper"].force).all()


def test_variable_hydro_wave_ramp_changes_excitation_before_state_switch():
    """The incident ramp applies to force and elevation, not the state clock."""
    hydro = BodyClass(str(SPHERE))
    hydro.bodyNumber = hydro.bodyTotal = 1
    hydro.readH5file()
    second_mass = 1.1 * 1000 * float(np.asarray(hydro.dispVol).item())
    settings = dict(height=2.5, switch_time=.16, end_time=.3)
    ramped = _variable_sphere(second_mass, ramp_time=.12, **settings)
    unramped = _variable_sphere(second_mass, ramp_time=0, **settings)

    time = ramped.time
    ramp = np.ones_like(time)
    early = time < .12
    ramp[early] = (1 - np.cos(np.pi * time[early] / .12)) / 2
    np.testing.assert_array_equal(ramped.time, unramped.time)
    np.testing.assert_allclose(ramped.wave_elevation,
                               unramped.wave_elevation * ramp,
                               rtol=0, atol=1e-13)
    for response in (ramped, unramped):
        state = dict(response.raw.extra_outputs)["body1_hydro_state"]
        np.testing.assert_array_equal(state[time < .16], 1)
        np.testing.assert_array_equal(state[time >= .16], 2)
        speed = response.ptos["damper"].velocity
        np.testing.assert_allclose(response.ptos["damper"].absorbed_power,
                                   200_000 * speed**2, rtol=0, atol=1e-9)

    force = dict(ramped.raw.extra_outputs)["body1_excitation_heave"]
    reference = dict(unramped.raw.extra_outputs)["body1_excitation_heave"]
    np.testing.assert_allclose(force, reference * ramp, rtol=1e-12, atol=1e-6)
    assert np.max(np.abs(ramped.bodies["sphere"].position[:, 2]
                         - unramped.bodies["sphere"].position[:, 2])) > 1e-7

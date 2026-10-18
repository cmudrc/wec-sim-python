"""Analytical checks for the PTO-Sim direct linear generator equations."""

from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from wecsim import RegularWave, WEC
from wecsim.components.directLinearGenerator import (
    DirectLinearGenerator, run_rm3_direct_linear_generator,
)


RM3 = (Path(__file__).parent / "test_objects" / "test_bodyclass" /
       "testData" / "hydroData" / "rm3.h5")


def test_constant_speed_flux_and_power_balance():
    generator = DirectLinearGenerator(
        stator_resistance=4.58, friction=-100, pole_pitch=0.072,
        magnet_flux=8, inductance=0.285, load_resistance=-117.6471,
    )
    speed = 0.4
    omega = np.pi * speed / generator.pole_pitch
    decay = ((generator.load_resistance - generator.stator_resistance)
             / generator.inductance)
    steady = np.linalg.solve(
        np.array([[decay, omega], [-omega, decay]]),
        np.array([decay * generator.magnet_flux, 0]),
    )
    # At constant speed, the two flux states rotate and decay toward a
    # closed-form equilibrium; angle advances at a constant electrical rate.
    duration = 0.05
    initial = generator.initial_state()
    rotation = np.array([
        [np.cos(omega * duration), np.sin(omega * duration)],
        [-np.sin(omega * duration), np.cos(omega * duration)],
    ])
    expected_flux = steady + np.exp(decay * duration) * rotation @ (
        initial[:2] - steady
    )
    solution = solve_ivp(
        lambda _t, state: generator.state_rate(speed, state),
        (0, duration), initial, rtol=1e-11, atol=1e-13,
    )
    np.testing.assert_allclose(solution.y[:2, -1], expected_flux,
                               rtol=1e-9, atol=1e-10)
    np.testing.assert_allclose(solution.y[2, -1], omega * duration,
                               rtol=1e-11, atol=1e-12)
    result = generator.signals(speed, solution.y[:, -1])
    assert result.force < 0
    assert result.absorbed_power > 0
    assert result.electrical_power > 0
    np.testing.assert_allclose(
        result.phase_current @ result.phase_current,
        result.current_d**2 + result.current_q**2,
        rtol=1e-13, atol=1e-13,
    )
    np.testing.assert_allclose(
        result.electrical_power,
        -generator.load_resistance
        * (result.current_d**2 + result.current_q**2),
        rtol=1e-13, atol=1e-10,
    )
    flux_rate = generator.state_rate(speed, solution.y[:, -1])
    magnetic_energy_rate = (
        result.current_d * flux_rate[0] + result.current_q * flux_rate[1]
    )
    stator_loss = generator.stator_resistance * (
        result.current_d**2 + result.current_q**2
    )
    np.testing.assert_allclose(
        -result.electromagnetic_force * speed,
        magnetic_energy_rate + stator_loss + result.electrical_power,
        rtol=1e-12, atol=1e-9,
    )


def test_unforced_rm3_generator_stays_at_equilibrium():
    generator = DirectLinearGenerator(
        stator_resistance=4.58, friction=-100, pole_pitch=0.072,
        magnet_flux=8, inductance=0.285, load_resistance=-117.6471,
    )
    response = run_rm3_direct_linear_generator(
        RM3, generator=generator, height=0, ramp_time=0,
        end_time=1, dt=0.0005, output_dt=0.01,
    )
    assert response.time.shape == (101,)
    np.testing.assert_allclose(response.body_heave,
                               np.broadcast_to(response.body_heave[0],
                                               response.body_heave.shape),
                               rtol=0, atol=1e-12)
    for values in (response.body_heave_velocity, response.pto_force,
                   response.absorbed_power, response.electrical_power,
                   response.phase_current, response.phase_voltage):
        np.testing.assert_allclose(values, 0, rtol=0, atol=1e-10)


def test_ramped_rm3_generator_agrees_with_public_wec():
    """Keep the standalone source-comparison solver paired with WEC.run."""
    settings = dict(
        stator_resistance=4.58, friction=-100, pole_pitch=.072,
        magnet_flux=8, inductance=.285, load_resistance=-117.6471,
    )
    source_solver = run_rm3_direct_linear_generator(
        RM3, generator=DirectLinearGenerator(**settings),
        height=2.5, period=8, ramp_time=100,
        dt=.0005, output_dt=.01, end_time=1,
    )

    wec = WEC("RM3 ramped direct generator")
    floating = wec.body("float", RM3)
    spar = wec.body("spar", RM3)
    wec.coordinate("float_heave", floating.move("heave"))
    wec.coordinate("spar_heave", spar.move("heave"))
    wec.pto(
        "PTO1", spar.at(0, 0, 0), floating.at(0, 0, 0),
        axis=(0, 0, 1),
        linear_generator=DirectLinearGenerator(**settings),
    )
    public = wec.run(RegularWave(2.5, 8), dt=.0005,
                     end_time=1, ramp_time=100)
    sample = slice(None, None, 20)  # 0.0005 s integration, 0.01 s output.
    np.testing.assert_allclose(source_solver.time, public.time[sample],
                               rtol=0, atol=1e-12)
    np.testing.assert_allclose(source_solver.wave_elevation,
                               public.wave_elevation[sample], rtol=0, atol=1e-12)
    for index, name in enumerate(("float", "spar")):
        np.testing.assert_allclose(
            source_solver.body_heave[:, index],
            public.bodies[name].position[sample, 2], rtol=0, atol=1e-12,
        )
        np.testing.assert_allclose(
            source_solver.body_heave_velocity[:, index],
            public.bodies[name].velocity[sample, 2], rtol=0, atol=1e-12,
        )
    pto = public.ptos["PTO1"]
    np.testing.assert_allclose(source_solver.pto_force, pto.force[sample],
                               rtol=0, atol=1e-10)
    np.testing.assert_allclose(source_solver.absorbed_power,
                               pto.absorbed_power[sample], rtol=0, atol=1e-10)

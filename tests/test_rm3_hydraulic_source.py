"""Pinned RM3 hydraulic PTO source signals and replayed incident sea.

This is an input/force diagnostic. It does not claim coupled hydraulic
trajectory parity with a Python PTO model.
"""

import os
from functools import cache
from pathlib import Path

import numpy as np
import pytest

from wecsim.components.electricGenerator import (
    DiscretePILoadController, EquivalentCircuitGenerator,
)
from wecsim.components.hydraulic import (
    CompressibleCylinder, ConstantEfficiencyHydraulicMotor,
    GasChargedAccumulator, RectifyingCheckValve,
)
from wecsim.environment.irregularWave import (
    pm_equal_energy_components, synthesize_irregular_response,
)


REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
RM3_H5 = os.environ.get("WEC_SIM_RM3_H5")
LIVE = bool(REFERENCE and RM3_H5)
FIXTURE = Path(__file__).parent / "fixtures/matlab/rm3_hydraulic"
HYDRO = Path(RM3_H5) if LIVE else FIXTURE / "rm3.h5"


@cache
def _load(name):
    if LIVE:
        values = np.loadtxt(Path(REFERENCE) / f"RM3_HYDRAULIC_SOURCE_{name}.csv",
                            delimiter=",")
    else:
        with np.load(FIXTURE / "source.npz", allow_pickle=False) as saved:
            values = saved[name]
    assert values.ndim == 2 and np.isfinite(values).all(), name
    return values


def _published_valve():
    area_max, area_min = .002, 1e-8
    pressure_max, pressure_min = 1.5e6, 0
    opening_gain = np.arctanh(
        (area_min - (area_max - area_min) / 2)
        * 2 / (area_max - area_min)
    ) / (pressure_min - (pressure_max + pressure_min) / 2)
    return RectifyingCheckValve(
        .61, area_max, area_min, pressure_max, pressure_min,
        850, 200, opening_gain,
    )


def test_hydraulic_source_wave_and_component_traces():
    components = _load("components")
    matlab_wave = _load("wave")
    assert components.shape == (500, 4)
    assert matlab_wave.shape == (40_001, 2)
    generated = pm_equal_energy_components(
        HYDRO, significant_height=2.5, peak_period=8,
        directions=(0,), spreading=(1,), phase=components[:, 3:4],
    )
    for actual, saved in ((generated.omega, components[:, 0]),
                          (generated.spectral_amplitude, components[:, 1]),
                          (generated.d_omega, components[:, 2])):
        assert np.max(np.abs(actual - saved)) < 1e-12
    replay = synthesize_irregular_response(
        HYDRO, generated, dt=.01, end_time=400, ramp_time=100,
    )
    np.testing.assert_allclose(replay.time, matlab_wave[:, 0], atol=1e-10,
                               rtol=0)
    assert np.max(np.abs(replay.elevation - matlab_wave[:, 1])) < 1e-10

    for name, columns in (("cylinder", 4), ("valve", 5),
                          ("high_accumulator", 3),
                          ("low_accumulator", 3), ("motor", 5),
                          ("generator", 5)):
        values = _load(name)
        assert values.shape == (40_001, columns), name
        np.testing.assert_allclose(values[:, 0], matlab_wave[:, 0],
                                   atol=1e-10, rtol=0)


def test_hydraulic_cylinder_pressure_and_force_balance():
    cylinder = _load("cylinder")
    valve = _load("valve")
    pto = _load("RM3_cHydraulic_PTO_pto1")
    assert pto.shape == (40_001, 25)
    np.testing.assert_allclose(pto[:, 0], cylinder[:, 0], atol=1e-10,
                               rtol=0)
    assert np.max(np.abs(cylinder[:, 2])) > 1e3
    model = CompressibleCylinder(.0378, .0378, 1.86e9, 70, 35)
    reconstructed = model.force(cylinder[:, 1], cylinder[:, 3])
    assert np.max(np.abs(reconstructed - cylinder[:, 2])) < 1e-4

    # The source piston selects heave from the PTO response bus. Its first
    # two valve outputs are chamber port flows, and pressure uses a discrete
    # integrator at the model time step.
    rate_a, rate_b = model.pressure_rates(
        pto[:, 3], pto[:, 9], valve[:, 1], valve[:, 2],
    )
    dt = np.diff(cylinder[:, 0])
    for pressure, rate in ((cylinder[:, 1], rate_a),
                           (cylinder[:, 3], rate_b)):
        predicted = pressure[:-1] + dt * rate[:-1]
        assert np.max(np.abs(predicted - pressure[1:])) < 3e-6


def test_hydraulic_rectifying_valve_port_flows():
    cylinder = _load("cylinder")
    valve_output = _load("valve")
    high = _load("high_accumulator")
    low = _load("low_accumulator")
    pto = _load("RM3_cHydraulic_PTO_pto1")
    model = _published_valve()
    flows = np.column_stack(model.flows(
        cylinder[:, 1], cylinder[:, 3], high[:, 1], low[:, 1],
    ))
    assert np.max(np.abs(valve_output[:, 1:])) > .01
    assert np.max(np.abs(flows - valve_output[:, 1:])) < 1e-12

    piston = CompressibleCylinder(.0378, .0378, 1.86e9, 70, 35)
    rate_a, rate_b = piston.pressure_rates(
        pto[:, 3], pto[:, 9], flows[:, 0], flows[:, 1],
    )
    dt = np.diff(cylinder[:, 0])
    for pressure, rate in ((cylinder[:, 1], rate_a),
                           (cylinder[:, 3], rate_b)):
        assert np.max(np.abs(pressure[:-1] + dt * rate[:-1]
                             - pressure[1:])) < 3e-6


@pytest.mark.parametrize(
    ("name", "precharge"),
    (("high_accumulator", 2784.7 * 6894.75),
     ("low_accumulator", 1392.4 * 6894.75)),
)
def test_hydraulic_accumulator_pressure_from_port_flow(name, precharge):
    trace = _load(name)
    cylinder = _load("cylinder")
    high = _load("high_accumulator")
    low = _load("low_accumulator")
    motor = _load("motor")
    flows = _published_valve().flows(
        cylinder[:, 1], cylinder[:, 3], high[:, 1], low[:, 1],
    )
    if name == "high_accumulator":
        inlet_flow = flows[2] - motor[:, 4]
    else:
        inlet_flow = flows[3] + motor[:, 4]
    model = GasChargedAccumulator(8.5, precharge)
    accepted_volume = np.r_[
        0, np.cumsum(np.diff(trace[:, 0]) * inlet_flow[:-1]),
    ]
    reconstructed = model.pressure(accepted_volume)
    assert np.max(np.abs(trace[:, 2])) > .01
    assert np.max(np.abs(inlet_flow - trace[:, 2])) < 1e-12
    assert np.max(np.abs(reconstructed - trace[:, 1])) < 1e-6


def test_hydraulic_motor_and_generator_drive_states():
    high = _load("high_accumulator")
    low = _load("low_accumulator")
    motor_output = _load("motor")
    generator_output = _load("generator")
    motor = ConstantEfficiencyHydraulicMotor(120e-6, .9, .85)
    generator = EquivalentCircuitGenerator(.8, .8, .8, .8, .8)
    controller = DiscretePILoadController(1000, .001, .001)

    pressure_drop = motor_output[:, 3]
    speed = generator_output[:, 2] * (2 * np.pi / 60)
    drive_torque = motor.torque(pressure_drop)
    motor_flow = motor.flow(speed)
    assert np.max(np.abs(high[:, 1] - low[:, 1] - pressure_drop)) < 2e-7
    assert np.max(np.abs(speed - motor_output[:, 1] * (2 * np.pi / 60))) < 1e-12
    assert np.max(np.abs(drive_torque - motor_output[:, 2])) < 1e-10
    assert np.max(np.abs(motor_flow - motor_output[:, 4])) < 1e-14

    current = generator_output[:, 3]
    dt = np.diff(generator_output[:, 0])
    integral = np.r_[
        0, np.cumsum(dt * controller.integral_rate(generator_output[:-1, 2])),
    ]
    voltage = controller.voltage(generator_output[:, 2], current, integral)
    assert np.max(np.abs(voltage - generator_output[:, 4])) < 1e-10
    torque_em = generator.electromagnetic_torque(current)
    assert np.max(np.abs(torque_em - generator_output[:, 1])) < 1e-10
    current_rate = generator.current_rate(speed, current, voltage)
    speed_rate = generator.speed_rate(speed, current, drive_torque)
    assert np.max(np.abs(current[:-1] + dt * current_rate[:-1]
                         - current[1:])) < 1e-10
    assert np.max(np.abs(speed[:-1] + dt * speed_rate[:-1]
                         - speed[1:])) < 1e-10

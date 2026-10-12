"""Independent 400 s body/PTO trajectory against the pinned RM3 hydraulic case."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim import (
    CompressibleCylinder, ConstantEfficiencyHydraulicMotor,
    DiscretePILoadController, EquivalentCircuitGenerator,
    GasChargedAccumulator, RectifiedHydraulicPTO, RectifyingCheckValve,
    run_rm3_rectified_hydraulic,
)
from wecsim.environment.irregularWave import pm_equal_energy_components


REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
RM3_H5 = os.environ.get("WEC_SIM_RM3_H5")
pytestmark = pytest.mark.skipif(
    not (REFERENCE and RM3_H5),
    reason="pinned RM3 hydraulic MATLAB output and HDF5 absent",
)


def _source(name):
    data = np.loadtxt(
        Path(REFERENCE) / f"RM3_HYDRAULIC_SOURCE_{name}.csv", delimiter=",",
    )
    assert data.ndim == 2 and np.isfinite(data).all(), name
    return data


def _published_pto():
    area_max, area_min = .002, 1e-8
    pressure_max, pressure_min = 1.5e6, 0.
    opening_gain = np.arctanh(
        (area_min - (area_max - area_min) / 2) * 2 / (area_max - area_min)
    ) / (pressure_min - (pressure_max + pressure_min) / 2)
    return RectifiedHydraulicPTO(
        cylinder=CompressibleCylinder(.0378, .0378, 1.86e9, 70, 35),
        valve=RectifyingCheckValve(
            .61, area_max, area_min, pressure_max, pressure_min,
            850, 200, opening_gain,
        ),
        high_accumulator=GasChargedAccumulator(8.5, 2784.7 * 6894.75),
        low_accumulator=GasChargedAccumulator(8.5, 1392.4 * 6894.75),
        motor=ConstantEfficiencyHydraulicMotor(120e-6, .9, .85),
        generator=EquivalentCircuitGenerator(.8, .8, .8, .8, .8),
        load_controller=DiscretePILoadController(1000, .001, .001),
        initial_pressure_a=2.1333e7,
        initial_pressure_b=2.1333e7,
    )


@pytest.fixture(scope="module")
def paired():
    components = _source("components")
    sea = pm_equal_energy_components(
        RM3_H5, significant_height=2.5, peak_period=8,
        directions=(0,), spreading=(1,), phase=components[:, 3:4],
    )
    response = run_rm3_rectified_hydraulic(RM3_H5, sea, _published_pto())
    source = {name: _source(name) for name in (
        "wave", "RM3_cHydraulic_PTO_body1", "RM3_cHydraulic_PTO_body2",
        "RM3_cHydraulic_PTO_pto1", "cylinder", "valve",
        "high_accumulator", "low_accumulator", "motor", "generator",
    )}
    assert response.time.shape == (40_001,)
    for name, values in source.items():
        assert values.shape[0] == len(response.time), name
        np.testing.assert_allclose(response.time, values[:, 0], rtol=0,
                                   atol=1e-9, err_msg=name)
    return response, source


def _bounded(actual, expected, bound, label):
    assert actual.shape == expected.shape, label
    assert np.isfinite(actual).all(), label
    error = float(np.max(np.abs(actual - expected)))
    assert error < bound, f"{label}: maximum {error:.6g} exceeds {bound}"


def test_published_two_heave_motion_and_incident_wave(paired):
    response, source = paired
    _bounded(response.wave_elevation, source["wave"][:, 1],
             1e-10, "incident wave")
    for index in (1, 2):
        body = source[f"RM3_cHydraulic_PTO_body{index}"]
        # The published hydraulic model's constrained bodies have only heave.
        assert np.max(np.abs(body[:, [1, 2, 4, 5, 6, 7, 8, 10, 11, 12]])) < 1e-10
        _bounded(response.body_heave[:, index - 1], body[:, 3],
                 .005, f"body {index} heave")
        _bounded(response.body_heave_velocity[:, index - 1], body[:, 9],
                 .004, f"body {index} heave velocity")
    pto = source["RM3_cHydraulic_PTO_pto1"]
    _bounded(response.pto_stroke, pto[:, 3], .006, "PTO stroke")
    _bounded(response.pto_velocity, pto[:, 9], .005, "PTO speed")


def test_coupled_hydraulic_and_electric_outputs(paired):
    response, source = paired
    cylinder = source["cylinder"]
    _bounded(response.pto_force, cylinder[:, 2], 20_000, "cylinder force")
    _bounded(response.cylinder_pressure[:, 0], cylinder[:, 1],
             400_000, "chamber A pressure")
    _bounded(response.cylinder_pressure[:, 1], cylinder[:, 3],
             400_000, "chamber B pressure")
    for index, name in enumerate(("high_accumulator", "low_accumulator")):
        _bounded(response.accumulator_pressure[:, index], source[name][:, 1],
                 2_000, f"{name} pressure")
    valve_error = response.valve_flow - source["valve"][:, 1:]
    _bounded(response.valve_flow, source["valve"][:, 1:], .05,
             "valve port flow")
    assert np.max(np.abs(np.cumsum(valve_error, axis=0) * .01)) < 6e-4
    motor = source["motor"]
    generator = source["generator"]
    _bounded(response.motor_torque, motor[:, 2], .04, "motor torque")
    _bounded(response.motor_flow, motor[:, 4], 1.5e-7, "motor flow")
    _bounded(response.shaft_speed_rpm, motor[:, 1], .06, "shaft speed")
    _bounded(response.generator_current, generator[:, 3], .03,
             "generator current")
    _bounded(response.generator_voltage, generator[:, 4], .03,
             "generator voltage")


def test_published_controller_is_source_specific(paired):
    response, _ = paired
    # This is the published PI load law, whose negative resistance is a
    # limitation of the source case rather than a passive-device assertion.
    assert response.load_resistance[0] > 0
    assert np.any(response.load_resistance < 0)

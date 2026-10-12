"""Independent 400 s pitch/PTO trajectories for both pinned OSWEC layouts."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim import (
    AdjustableRodCrank, CompressibleCylinder,
    ConstantEfficiencyHydraulicMotor, DiscretePILoadController,
    EquivalentCircuitGenerator, FixedRodCrank, GasChargedAccumulator,
    RectifiedHydraulicPTO, RectifyingCheckValve,
    run_oswec_rectified_hydraulic,
)
from wecsim.environment.irregularWave import pm_equal_energy_components


REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
OSWEC_H5 = os.environ.get("WEC_SIM_OSWEC_H5")
pytestmark = pytest.mark.skipif(
    not (REFERENCE and OSWEC_H5),
    reason="pinned OSWEC hydraulic MATLAB output and HDF5 absent",
)


def _source(case, name):
    values = np.loadtxt(
        Path(REFERENCE) / f"OSWEC_HYDRAULIC_SOURCE_{case}_{name}.csv",
        delimiter=",",
    )
    assert values.ndim == 2 and np.isfinite(values).all(), name
    return values


def _published_pto(fixed_crank):
    maximum, minimum = .002, 1e-8
    opening = np.arctanh(
        (minimum - (maximum - minimum) / 2) * 2 / (maximum - minimum)
    ) / (-1.5e6 / 2)
    return RectifiedHydraulicPTO(
        cylinder=CompressibleCylinder(.0378, .0378, 1.86e9, 70, 35),
        valve=RectifyingCheckValve(
            .61, maximum, minimum, 1.5e6, 0, 850, 200, opening,
        ),
        high_accumulator=GasChargedAccumulator(8.5, 2784.7 * 6894.75),
        low_accumulator=GasChargedAccumulator(8.5, 1392.4 * 6894.75),
        motor=ConstantEfficiencyHydraulicMotor(120e-6, .9, .85),
        generator=EquivalentCircuitGenerator(
            1, .1, 1, 2 if fixed_crank else 20, .1,
        ),
        load_controller=DiscretePILoadController(
            3000 if fixed_crank else 1500, .001, .001,
        ),
        initial_pressure_a=1.4e7,
        initial_pressure_b=1.4e7,
    )


@pytest.fixture(scope="module", params=(
    ("OSWEC_Hydraulic_PTO", False),
    ("OSWEC_Hydraulic_Crank_PTO", True),
))
def paired(request):
    case, fixed = request.param
    components = _source(case, "components")
    sea = pm_equal_energy_components(
        OSWEC_H5, significant_height=2.5, peak_period=8,
        directions=(0,), spreading=(1,), phase=components[:, 3:4],
    )
    linkage = FixedRodCrank(3, 1.3, 5) if fixed else AdjustableRodCrank(3, 1.3, 5)
    response = run_oswec_rectified_hydraulic(
        OSWEC_H5, sea, _published_pto(fixed), linkage,
    )
    source = {name: _source(case, name) for name in (
        "wave", "body1", "body2", "pto1", "crank", "cylinder",
        "valve", "high_accumulator", "low_accumulator", "motor",
        "generator",
    )}
    assert response.time.shape == (40_001,)
    for name, values in source.items():
        np.testing.assert_allclose(response.time, values[:, 0], rtol=0,
                                   atol=1e-9, err_msg=name)
    return case, response, source, linkage


def _bound(actual, expected, limit, label):
    assert actual.shape == expected.shape, label
    assert np.isfinite(actual).all(), label
    error = float(np.max(np.abs(actual - expected)))
    assert error < limit, f"{label}: maximum {error:.6g} exceeds {limit}"


def test_published_flap_and_crank_motion(paired):
    case, response, source, linkage = paired
    body = source["body1"]
    _bound(response.wave_elevation, source["wave"][:, 1],
           1e-10, f"{case} wave")
    _bound(response.excitation_force, body[:, 19:25],
           1e-6, f"{case} excitation")
    _bound(response.pitch, body[:, 5], .002, f"{case} pitch")
    _bound(response.pitch_velocity, body[:, 11], .001,
           f"{case} pitch speed")
    _bound(response.body_position, body[:, 1:4], .008,
           f"{case} flap center")
    _bound(response.body_velocity, body[:, 7:10], .005,
           f"{case} flap center speed")
    base = source["body2"]
    assert np.max(np.abs(base[:, 1:7] - base[0, 1:7])) < 1e-10
    assert np.max(np.abs(base[:, 7:13])) < 1e-10
    _bound(response.crank_torque, source["crank"][:, 1], 75_000,
           f"{case} crank torque")
    source_stroke, source_jacobian = linkage.stroke_and_jacobian(
        source["crank"][:, 2],
    )
    _bound(response.crank_stroke, source_stroke, .008,
           f"{case} crank stroke")
    _bound(response.crank_speed, source_jacobian * source["crank"][:, 3],
           .005, f"{case} crank speed")


def test_coupled_hydraulic_and_electrical_outputs(paired):
    case, response, source, _ = paired
    cylinder = source["cylinder"]
    _bound(response.cylinder_force, cylinder[:, 2], 25_000,
           f"{case} cylinder force")
    _bound(response.cylinder_pressure[:, 0], cylinder[:, 1], 500_000,
           f"{case} chamber A")
    _bound(response.cylinder_pressure[:, 1], cylinder[:, 3], 500_000,
           f"{case} chamber B")
    for index, name in enumerate(("high_accumulator", "low_accumulator")):
        _bound(response.accumulator_pressure[:, index], source[name][:, 1],
               7_000, f"{case} {name}")
    flow_error = response.valve_flow - source["valve"][:, 1:]
    _bound(response.valve_flow, source["valve"][:, 1:], .065,
           f"{case} valve flow")
    assert np.max(np.abs(np.cumsum(flow_error, axis=0) * .01)) < .002
    _bound(response.shaft_speed_rpm, source["motor"][:, 1], .08,
           f"{case} shaft speed")
    _bound(response.motor_torque, source["motor"][:, 2], .15,
           f"{case} motor torque")
    _bound(response.motor_flow, source["motor"][:, 4], 2e-7,
           f"{case} motor flow")
    _bound(response.generator_current, source["generator"][:, 3], .2,
           f"{case} current")
    _bound(response.generator_voltage, source["generator"][:, 4], .25,
           f"{case} voltage")

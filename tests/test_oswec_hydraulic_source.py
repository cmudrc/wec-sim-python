"""Source-driven PTO checks for both pinned OSWEC hydraulic applications."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim.components.crank import AdjustableRodCrank, FixedRodCrank
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
OSWEC_H5 = os.environ.get("WEC_SIM_OSWEC_H5")
pytestmark = pytest.mark.skipif(
    not (REFERENCE and OSWEC_H5),
    reason="pinned OSWEC hydraulic MATLAB output and HDF5 absent",
)
CASES = (
    ("OSWEC_Hydraulic_PTO", AdjustableRodCrank(3, 1.3, 5), 20, 1500, 1),
    ("OSWEC_Hydraulic_Crank_PTO", FixedRodCrank(3, 1.3, 5), 2, 3000, -1),
)


def _load(case, name):
    values = np.loadtxt(
        Path(REFERENCE) / f"OSWEC_HYDRAULIC_SOURCE_{case}_{name}.csv",
        delimiter=",",
    )
    assert values.ndim == 2 and np.isfinite(values).all(), name
    return values


def _valve():
    maximum, minimum = .002, 1e-8
    opening = np.arctanh(
        (minimum - (maximum - minimum) / 2) * 2 / (maximum - minimum)
    ) / (-1.5e6 / 2)
    return RectifyingCheckValve(
        .61, maximum, minimum, 1.5e6, 0, 850, 200, opening,
    )


@pytest.mark.parametrize("case,linkage,inertia,reference_rpm,torque_sign", CASES)
def test_oswec_hydraulic_wave_and_crank(case, linkage, inertia,
                                        reference_rpm, torque_sign):
    components = _load(case, "components")
    wave = _load(case, "wave")
    body = _load(case, "body1")
    crank = _load(case, "crank")
    cylinder = _load(case, "cylinder")
    assert components.shape == (500, 4)
    assert wave.shape == (40_001, 2)
    assert body.shape == (40_001, 25)
    assert crank.shape == (40_001, 4)
    sea = pm_equal_energy_components(
        OSWEC_H5, significant_height=2.5, peak_period=8,
        directions=(0,), spreading=(1,), phase=components[:, 3:4],
    )
    for actual, saved in ((sea.omega, components[:, 0]),
                          (sea.spectral_amplitude, components[:, 1]),
                          (sea.d_omega, components[:, 2])):
        assert np.max(np.abs(actual - saved)) < 1e-12
    excitation = synthesize_irregular_response(
        OSWEC_H5, sea, dt=.01, end_time=400, ramp_time=100,
    )
    assert np.max(np.abs(excitation.elevation - wave[:, 1])) < 1e-10
    assert np.max(np.abs(excitation.excitation_force[:, :6]
                         - body[:, 19:25])) < 1e-6
    assert np.max(np.abs(crank[:, 2] - body[:, 5])) < 1e-12
    assert np.max(np.abs(crank[:, 3] - body[:, 11])) < 1e-12
    _, jacobian = linkage.stroke_and_jacobian(crank[:, 2])
    calculated_torque = linkage.torque(cylinder[:, 2], crank[:, 2])
    assert np.max(np.abs(calculated_torque
                         - torque_sign * cylinder[:, 2] * jacobian)) < 1e-9
    assert np.max(np.abs(crank[:, 1] - calculated_torque)) < 1e-5


@pytest.mark.parametrize("case,linkage,inertia,reference_rpm,torque_sign", CASES)
def test_oswec_hydraulic_network_steps(case, linkage, inertia,
                                       reference_rpm, torque_sign):
    crank = _load(case, "crank")
    cylinder = _load(case, "cylinder")
    valve = _load(case, "valve")
    high = _load(case, "high_accumulator")
    low = _load(case, "low_accumulator")
    motor = _load(case, "motor")
    generator = _load(case, "generator")
    assert cylinder.shape == (40_001, 4)
    assert valve.shape == (40_001, 5)
    assert motor.shape == generator.shape == (40_001, 5)
    dt = np.diff(cylinder[:, 0])
    stroke, jacobian = linkage.stroke_and_jacobian(crank[:, 2])
    speed = jacobian * crank[:, 3]
    piston = CompressibleCylinder(.0378, .0378, 1.86e9, 70, 35)
    assert np.max(np.abs(piston.force(cylinder[:, 1], cylinder[:, 3])
                         - cylinder[:, 2])) < 1e-5
    flows = np.column_stack(_valve().flows(
        cylinder[:, 1], cylinder[:, 3], high[:, 1], low[:, 1],
    ))
    assert np.max(np.abs(flows - valve[:, 1:])) < 1e-12
    rate_a, rate_b = piston.pressure_rates(
        stroke, speed, flows[:, 0], flows[:, 1],
    )
    for pressure, rate in ((cylinder[:, 1], rate_a),
                           (cylinder[:, 3], rate_b)):
        assert np.max(np.abs(pressure[:-1] + dt * rate[:-1]
                             - pressure[1:])) < 3e-6

    drive = ConstantEfficiencyHydraulicMotor(120e-6, .9, .85)
    omega = motor[:, 1] * 2 * np.pi / 60
    assert np.max(np.abs(drive.torque(high[:, 1] - low[:, 1])
                         - motor[:, 2])) < 1e-10
    assert np.max(np.abs(drive.flow(omega) - motor[:, 4])) < 1e-14
    for accumulator, trace, inlet in (
        (GasChargedAccumulator(8.5, 2784.7 * 6894.75), high,
         flows[:, 2] - motor[:, 4]),
        (GasChargedAccumulator(8.5, 1392.4 * 6894.75), low,
         flows[:, 3] + motor[:, 4]),
    ):
        volume = np.r_[0, np.cumsum(dt * inlet[:-1])]
        assert np.max(np.abs(accumulator.pressure(volume)
                             - trace[:, 1])) < 1e-6

    electric = EquivalentCircuitGenerator(1, .1, 1, inertia, .1)
    controller = DiscretePILoadController(reference_rpm, .001, .001)
    integral = np.r_[0, np.cumsum(dt * controller.integral_rate(
        generator[:-1, 2]))]
    assert np.max(np.abs(controller.voltage(
        generator[:, 2], generator[:, 3], integral
    ) - generator[:, 4])) < 1e-8
    assert np.max(np.abs(generator[:-1, 3] + dt * electric.current_rate(
        omega[:-1], generator[:-1, 3], generator[:-1, 4],
    ) - generator[1:, 3])) < 1e-9
    assert np.max(np.abs(omega[:-1] + dt * electric.speed_rate(
        omega[:-1], generator[:-1, 3], motor[:-1, 2],
    ) - omega[1:])) < 1e-9

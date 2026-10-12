"""Independent high-pressure network against pinned 300 s Desalination."""

import os
from pathlib import Path

import numpy as np
import pytest
from scipy.ndimage import uniform_filter1d

from wecsim import (
    DynamicPressureReliefValve, FourValveRectifiedCylinder,
    GasChargedAccumulator, IdealDoubleActingCylinder,
    MorisonElement, PitchRodLinkage, ReverseOsmosisHydraulicNetwork,
    ReverseOsmosisMembrane, run_oswec_desalination,
)
from wecsim.environment.irregularWave import pm_equal_energy_components


REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
pytestmark = pytest.mark.skipif(
    not REFERENCE, reason="pinned OSWEC Desalination output absent",
)


def _source(name):
    return np.loadtxt(
        Path(REFERENCE) / f"OSWEC_DESALINATION_SOURCE_{name}.csv",
        delimiter=",",
    )


def _network():
    return ReverseOsmosisHydraulicNetwork(
        accumulator=GasChargedAccumulator(
            initial_gas_volume=4, precharge_pressure=3e6,
            atmospheric_pressure=101_325, dead_gas_volume=4e-5,
            hard_stop_stiffness=1e7, hard_stop_damping=1e7,
        ),
        membrane=ReverseOsmosisMembrane(
            resistance=.6023e8, osmotic_pressure=3e6,
            regulation_range=5e4, max_area=.3, leakage_area=1e-12,
            discharge_coefficient=.7, fluid_density=850,
        ),
        relief=DynamicPressureReliefValve(
            set_pressure=5.6e6, regulation_range=1e5,
            max_area=.0267, leakage_area=1e-12,
            discharge_coefficient=.7, fluid_density=850,
            opening_time_constant=.1,
        ),
        cylinder_area=.26,
        pump_to_motor_displacement_ratio=.95,
        pump_outlet_resistance=.6023e8 / 22,
    )


def test_published_relief_and_pressure_recovery_laws():
    network = _network()
    hydraulic, relief = _source("simout1"), _source("relief")
    pressure, time = hydraulic[:, 6], hydraulic[:, 0]
    area = np.empty_like(pressure)
    area[0] = relief[0, 2]
    for index in range(1, len(pressure)):
        area[index] = network.relief.next_area(
            pressure[index], area[index - 1], time[index] - time[index - 1],
        )
    np.testing.assert_allclose(area, relief[:, 2], rtol=0, atol=1e-12)
    actual_relief = np.array([
        network.relief.flow(p, a) for p, a in zip(pressure, area)
    ])
    np.testing.assert_allclose(actual_relief, relief[:, 3], rtol=0, atol=1e-10)
    ratio = network.pump_to_motor_displacement_ratio
    brine = pressure * (1 - ratio) / (
        ratio * ratio * network.pump_outlet_resistance
    )
    np.testing.assert_allclose(brine, hydraulic[:, 3], rtol=0, atol=1e-12)
    np.testing.assert_allclose(ratio * brine, hydraulic[:, 4],
                               rtol=0, atol=1e-12)


def test_published_network_from_source_feed_flow():
    hydraulic = _source("simout1")
    network = _network()
    state = network.initial_state(pressure=hydraulic[0, 6])
    predicted = np.empty(6001)
    predicted[0] = state.pressure
    for index in range(1, len(predicted)):
        state = network.step_from_feed_flow(
            hydraulic[index, 1], state,
            hydraulic[index, 0] - hydraulic[index - 1, 0],
        )
        predicted[index] = state.pressure
    np.testing.assert_allclose(predicted, hydraulic[:len(predicted), 6],
                               rtol=0, atol=1e-2)


@pytest.fixture(scope="module")
def rod_driven_trace():
    hydraulic, cylinder, relief = (
        _source("simout1"), _source("cylinder"), _source("relief")
    )
    assert hydraulic.shape == (30_001, 8)
    assert cylinder.shape == (30_001, 8)
    network = _network()
    state = network.initial_state()
    predicted = np.empty((len(hydraulic), 6))
    predicted[0] = 0
    for index in range(1, len(hydraulic)):
        state = network.step(
            cylinder[index, 2], state,
            hydraulic[index, 0] - hydraulic[index - 1, 0],
        )
        predicted[index] = (
            state.pressure, state.permeate_flow, state.brine_flow,
            state.recovered_feed_flow, state.accumulator_flow,
            state.relief_flow,
        )
    return hydraulic, cylinder, relief, predicted


def test_published_network_from_rod_speed_alone(rod_driven_trace):
    hydraulic, _, relief, predicted = rod_driven_trace
    references = np.column_stack((
        hydraulic[:, 6], hydraulic[:, 2], hydraulic[:, 3],
        hydraulic[:, 4], hydraulic[:, 7], relief[:, 3],
    ))
    limits = [500, 1e-5, 1e-5, 1e-5, 3e-3, 3e-3]
    assert np.all(np.max(np.abs(predicted - references), axis=0) < limits)
    assert np.max(predicted[:, 0]) > 5e6


def test_published_chamber_force_without_matching_pressure_chatter(
        rod_driven_trace):
    hydraulic, cylinder, _, predicted = rod_driven_trace
    valves = FourValveRectifiedCylinder(
        cylinder=IdealDoubleActingCylinder(.26, .26),
        max_area=.05, leakage_area=1e-8,
        discharge_coefficient=.7, fluid_density=850,
    )
    force = np.array([
        valves.rod_force(speed, pressure)
        for speed, pressure in zip(cylinder[:, 2], predicted[:, 0])
    ])
    source_force = cylinder[:, 7]
    # The pinned legacy solver alternates chamber pressure every 0.01 s;
    # smoothing both signals over 0.1 s tests their resolved force. The raw
    # mismatch remains explicit so this cannot be read as pointwise parity.
    raw_rms = np.sqrt(np.mean((force - source_force) ** 2))
    smooth_error = (
        uniform_filter1d(force, size=10, mode="nearest")
        - uniform_filter1d(source_force, size=10, mode="nearest")
    )
    assert raw_rms > 1e6
    assert np.sqrt(np.mean(smooth_error ** 2)) < 1e5
    source_work = np.trapezoid(source_force * cylinder[:, 2], hydraulic[:, 0])
    python_work = np.trapezoid(force * cylinder[:, 2], hydraulic[:, 0])
    assert abs(python_work - source_work) / abs(source_work) < .03

    # The valve log confirms why a raw-pressure gate would be misleading:
    # signed flow and logged pressure drop oppose each other in many samples.
    valve = _source("valve1")
    np.testing.assert_array_equal(valve[:, 0], hydraulic[:, 0])
    np.testing.assert_array_equal(valve[:, 1] > 0, cylinder[:, 2] > 0)
    assert np.mean(valve[:, 2] * valve[:, 3] < -1e-6) > .2


def test_published_fully_coupled_flap_and_desalination():
    hydro_file = os.environ.get("WEC_SIM_OSWEC_H5")
    if not hydro_file:
        pytest.skip("pinned OSWEC hydrodynamics absent")
    source_components = _source("components")
    components = pm_equal_energy_components(
        hydro_file, significant_height=2.64, peak_period=9.86,
        directions=(0,), spreading=(1,), count=250,
        phase=source_components[:, 3:4],
    )
    linkage = PitchRodLinkage(
        anchor=(5.6021271782, -8.7), hinge=(0, -8.9),
        body_center=(0, -3.9), body_point=(.9, -3.1),
    )
    elements = [MorisonElement(
        point=(0, 0, z), drag_coefficient=(1, 1, 1),
        added_mass_coefficient=(0, 0, 0), area=(32.4, 0, 32.4),
        volume=0,
    ) for z in (-3, -1.2, .6, 2.4, 4.2)]
    response = run_oswec_desalination(
        hydro_file, components, _network(),
        FourValveRectifiedCylinder(
            cylinder=IdealDoubleActingCylinder(.26, .26),
            max_area=.05, leakage_area=1e-8,
            discharge_coefficient=.7, fluid_density=850,
        ),
        linkage, elements,
    )
    body, cylinder, hydraulic, wave = (
        _source("Desalination_body1"), _source("cylinder"),
        _source("simout1"), _source("wave"),
    )
    assert response.time.shape == (30_001,)
    np.testing.assert_allclose(response.time, body[:, 0], rtol=0, atol=1e-12)
    np.testing.assert_allclose(response.wave_elevation, wave[:, 1],
                               rtol=0, atol=2e-12)
    np.testing.assert_allclose(response.excitation_force, body[:, 19:25],
                               rtol=0, atol=1e-6)
    assert np.max(np.abs(response.pitch - body[:, 5])) < .02
    assert np.max(np.abs(response.pitch_velocity - body[:, 11])) < .01
    assert np.max(np.abs(response.body_position - body[:, 1:4])) < .1
    assert np.max(np.abs(response.body_velocity - body[:, 7:10])) < .05
    assert np.max(np.abs(response.rod_speed - cylinder[:, 2])) < .02
    assert np.max(np.abs(response.high_pressure - hydraulic[:, 6])) < 40_000
    stroke, _ = linkage.stroke_and_jacobian(response.pitch)
    source_stroke, _ = linkage.stroke_and_jacobian(body[:, 5])
    assert np.max(np.abs(stroke - source_stroke)) < .05

    # A 0.1 s mean removes the legacy source's alternating chamber-pressure
    # artifact; total work remains a separate energy check.
    force_error = (
        uniform_filter1d(response.rod_force, size=10, mode="nearest")
        - uniform_filter1d(cylinder[:, 7], size=10, mode="nearest")
    )
    assert np.sqrt(np.mean(force_error ** 2)) < 150_000
    source_work = np.trapezoid(cylinder[:, 7] * cylinder[:, 2], body[:, 0])
    python_work = np.trapezoid(
        response.rod_force * response.rod_speed, response.time,
    )
    assert abs(python_work - source_work) / abs(source_work) < .02

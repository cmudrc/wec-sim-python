"""Short passive reverse-osmosis checks with bundled OSWEC hydrodynamics."""

from pathlib import Path

import numpy as np
import pytest

from wecsim.components.crank import PitchRodLinkage
from wecsim.components.desalination import (
    DynamicPressureReliefValve, FourValveRectifiedCylinder,
    ReverseOsmosisHydraulicNetwork, ReverseOsmosisMembrane,
)
from wecsim.components.hydraulic import (
    GasChargedAccumulator, IdealDoubleActingCylinder,
)
from wecsim.components.oswecDesalination import run_oswec_desalination
from wecsim.environment.irregularWave import pm_equal_energy_components
from wecsim.hydrodynamics.morison import MorisonElement


OSWEC = Path(__file__).resolve().parents[1] / "examples/data/oswec.h5"


def _hydraulics():
    network = ReverseOsmosisHydraulicNetwork(
        accumulator=GasChargedAccumulator(
            10, 1e6, hard_stop_stiffness=1e6,
        ),
        membrane=ReverseOsmosisMembrane(1e9, 1e5, 1e5, 1e-4,
                                        1e-8, .6, 1000),
        relief=DynamicPressureReliefValve(
            2e6, 1e5, 1e-4, 1e-8, .6, 1000, .1,
        ),
        cylinder_area=.02, pump_to_motor_displacement_ratio=.5,
        pump_outlet_resistance=1e10,
    )
    valves = FourValveRectifiedCylinder(
        IdealDoubleActingCylinder(.02, .02), 1e-4, 1e-8, .6, 1000,
    )
    return network, valves


def test_desalination_network_conserves_instantaneous_flow():
    network, valves = _hydraulics()
    initial = network.initial_state()
    assert abs(initial.pressure) < 1e-12
    step = network.step(.1, initial, .01)
    assert step.pressure >= 0
    assert step.permeate_flow >= 0
    np.testing.assert_allclose(
        step.accumulator_flow,
        step.cylinder_feed_flow + step.recovered_feed_flow
        - step.brine_flow - step.permeate_flow - step.relief_flow,
        rtol=0, atol=1e-12,
    )
    assert step.liquid_volume > initial.liquid_volume
    forward = valves.chamber_pressures(.1, step.pressure)
    reverse = valves.chamber_pressures(-.1, step.pressure)
    np.testing.assert_allclose(forward, reverse[::-1], rtol=0, atol=1e-7)
    np.testing.assert_allclose(valves.rod_force(.1, step.pressure),
                               .02 * (forward[0] - forward[1]))


def test_relief_valve_keeps_finite_opening_and_signed_flow():
    valve = _hydraulics()[0].relief
    area = valve.next_area(2.2e6, valve.leakage_area, .01)
    assert valve.leakage_area < area < valve.max_area
    assert valve.flow(0, area) == 0
    assert valve.flow(1e5, area) > 0
    assert valve.flow(-1e5, area) < 0

    for pressure, previous_area, dt in (
        (np.nan, area, .01), (2e6, np.inf, .01),
        (2e6, -area, .01), (2e6, area, 0),
    ):
        with pytest.raises(ValueError, match="relief opening"):
            valve.next_area(pressure, previous_area, dt)
    for pressure_drop, flow_area in ((np.nan, area), (1e5, -area)):
        with pytest.raises(ValueError, match="relief flow"):
            valve.flow(pressure_drop, flow_area)


def test_short_desalination_run_maps_rod_speed_and_force():
    network, valves = _hydraulics()
    linkage = PitchRodLinkage(
        anchor=(5.6021271782, -8.7), hinge=(0, -8.9),
        body_center=(0, -3.9), body_point=(.9, -3.1),
    )
    sea = pm_equal_energy_components(
        OSWEC, significant_height=.2, peak_period=8,
        directions=[0], spreading=[1], count=8, seed=2,
    )
    # A zero-drag element exercises the coupled entry point while keeping
    # the wave and hydraulic force checks independent of Morison tuning.
    element = MorisonElement((0, 0, 0), (0, 0, 0), (0, 0, 0),
                             (0, 0, 0), 0)
    result = run_oswec_desalination(
        OSWEC, sea, network, valves, linkage, [element],
        dt=.01, end_time=.05, ramp_time=0, radiation_memory=.02,
    )
    assert result.pitch.shape == (6,)
    assert np.isfinite(result.pitch).all()
    assert np.all(result.high_pressure >= 0)
    _, jacobian = linkage.stroke_and_jacobian(result.pitch)
    np.testing.assert_allclose(result.rod_speed,
                               -jacobian * result.pitch_velocity)
    np.testing.assert_allclose(result.rod_force,
                               valves.cylinder.force(*result.chamber_pressure.T))

"""Short coupled hydraulic runs with bundled RM3 and OSWEC hydrodynamics."""

from pathlib import Path

import numpy as np

from wecsim.components.crank import FixedRodCrank
from wecsim.components.electricGenerator import (
    DiscretePILoadController, EquivalentCircuitGenerator,
)
from wecsim.components.hydraulic import (
    CompressibleCylinder, ConstantEfficiencyHydraulicMotor,
    GasChargedAccumulator, RectifiedHydraulicPTO, RectifyingCheckValve,
)
from wecsim.components.oswecHydraulic import run_oswec_rectified_hydraulic
from wecsim.components.rm3Hydraulic import run_rm3_rectified_hydraulic
from wecsim.environment.irregularWave import pm_equal_energy_components


DATA = Path(__file__).resolve().parents[1] / "examples/data"


def _pto():
    return RectifiedHydraulicPTO(
        cylinder=CompressibleCylinder(.02, .02, 1e5, 10, 5),
        valve=RectifyingCheckValve(.6, 1e-6, 1e-8, 1e6, 0,
                                   850, 1e-5, 1e-5),
        high_accumulator=GasChargedAccumulator(10, 1e6),
        low_accumulator=GasChargedAccumulator(10, 5e5),
        motor=ConstantEfficiencyHydraulicMotor(1e-5, .9, .9),
        generator=EquivalentCircuitGenerator(1, 1, 1, 100, .1),
        load_controller=DiscretePILoadController(100, .001, .001),
        initial_pressure_a=7e5, initial_pressure_b=7e5,
    )


def _sea(path):
    return pm_equal_energy_components(
        path, significant_height=.2, peak_period=8,
        directions=[0], spreading=[1], count=8, seed=3,
    )


def test_rm3_hydraulic_runner_tracks_relative_stroke_and_pressure_force():
    hydro = DATA / "rm3.h5"
    pto = _pto()
    result = run_rm3_rectified_hydraulic(
        hydro, _sea(hydro), pto, dt=.01, end_time=.1,
        ramp_time=0, radiation_memory=.05,
    )
    assert result.body_heave.shape == (11, 2)
    assert np.isfinite(result.body_heave).all()
    np.testing.assert_allclose(result.pto_velocity,
                               result.body_heave_velocity[:, 0]
                               - result.body_heave_velocity[:, 1])
    np.testing.assert_allclose(result.pto_force,
                               pto.cylinder.force(result.cylinder_pressure[:, 0],
                                                  result.cylinder_pressure[:, 1]))
    np.testing.assert_allclose(result.generator_voltage,
                               result.load_resistance * result.generator_current)
    np.testing.assert_allclose(result.valve_flow.sum(axis=1), 0, atol=1e-12)


def test_oswec_hydraulic_runner_maps_cylinder_force_through_crank():
    hydro = DATA / "oswec.h5"
    pto = _pto()
    linkage = FixedRodCrank(3, 1.3, 5)
    result = run_oswec_rectified_hydraulic(
        hydro, _sea(hydro), pto, linkage,
        dt=.01, end_time=.1, ramp_time=0, radiation_memory=.05,
    )
    assert result.pitch.shape == (11,)
    assert np.isfinite(result.pitch).all()
    stroke, jacobian = linkage.stroke_and_jacobian(result.pitch)
    np.testing.assert_allclose(result.crank_stroke, stroke)
    np.testing.assert_allclose(result.crank_speed, jacobian * result.pitch_velocity)
    np.testing.assert_allclose(result.crank_torque,
                               -result.cylinder_force * jacobian)
    np.testing.assert_allclose(result.cylinder_force,
                               pto.cylinder.force(result.cylinder_pressure[:, 0],
                                                  result.cylinder_pressure[:, 1]))

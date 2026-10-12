"""Configure a published OSWEC hydraulic PTO entirely in Python.

Run ``python -m examples.oswec_hydraulic path/to/oswec.h5``. Choose
``--layout fixed`` for the fixed-rod crank; adjustable is the default.
"""

import argparse
from pathlib import Path

from wecsim import (
    AdjustableRodCrank, CompressibleCylinder,
    ConstantEfficiencyHydraulicMotor, DiscretePILoadController,
    EquivalentCircuitGenerator, FixedRodCrank, GasChargedAccumulator,
    RectifiedHydraulicPTO, RectifyingCheckValve,
    run_oswec_rectified_hydraulic,
)
from wecsim.environment.irregularWave import pm_equal_energy_components


def build_pto(*, fixed_crank: bool = False) -> RectifiedHydraulicPTO:
    return RectifiedHydraulicPTO(
        cylinder=CompressibleCylinder(
            area_a=.0378, area_b=.0378, bulk_modulus=1.86e9,
            stroke=70, offset=35,
        ),
        valve=RectifyingCheckValve(
            discharge_coefficient=.61, area_max=.002, area_min=1e-8,
            pressure_max=1.5e6, pressure_min=0, density=850,
            switch_gain=200, opening_gain=8.137375096991404e-6,
        ),
        high_accumulator=GasChargedAccumulator(8.5, 2784.7 * 6894.75),
        low_accumulator=GasChargedAccumulator(8.5, 1392.4 * 6894.75),
        motor=ConstantEfficiencyHydraulicMotor(120e-6, .9, .85),
        generator=EquivalentCircuitGenerator(
            armature_resistance=1, armature_inductance=.1,
            torque_constant=1, rotor_inertia=2 if fixed_crank else 20,
            shaft_damping=.1,
        ),
        load_controller=DiscretePILoadController(
            reference_rpm=3000 if fixed_crank else 1500,
            proportional_gain=.001, integral_gain=.001,
        ),
        initial_pressure_a=1.4e7,
        initial_pressure_b=1.4e7,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("hydro_file", type=Path)
    parser.add_argument("--layout", choices=("adjustable", "fixed"),
                        default="adjustable")
    args = parser.parse_args()
    fixed = args.layout == "fixed"
    linkage = (FixedRodCrank(3, 1.3, 5) if fixed
               else AdjustableRodCrank(3, 1.3, 5))
    sea = pm_equal_energy_components(
        args.hydro_file, significant_height=2.5, peak_period=8,
        directions=(0,), spreading=(1,), seed=1,
    )
    result = run_oswec_rectified_hydraulic(
        args.hydro_file, sea, build_pto(fixed_crank=fixed), linkage,
        dt=.01, end_time=20, ramp_time=5,
    )
    print(f"Samples: {len(result.time)}")
    print(f"Peak flap pitch: {abs(result.pitch).max():.4f} rad")
    print(f"Peak crank torque: {abs(result.crank_torque).max():.0f} N m")

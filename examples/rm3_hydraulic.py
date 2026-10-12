"""Configure and run the published two-heave RM3 hydraulic PTO in Python.

Run from the repository root with ``python -m examples.rm3_hydraulic``.
The seed gives a reproducible Python sea; MATLAB's phase realization differs.
"""

from pathlib import Path

from wecsim import (
    CompressibleCylinder, ConstantEfficiencyHydraulicMotor,
    DiscretePILoadController, EquivalentCircuitGenerator,
    GasChargedAccumulator, RectifiedHydraulicPTO, RectifyingCheckValve,
    run_rm3_rectified_hydraulic,
)
from wecsim.environment.irregularWave import pm_equal_energy_components


HYDRO = Path(__file__).resolve().parent / "data/rm3.h5"


def build_pto() -> RectifiedHydraulicPTO:
    # All settings are ordinary Python values and can be changed before run.
    return RectifiedHydraulicPTO(
        cylinder=CompressibleCylinder(
            area_a=.0378, area_b=.0378, bulk_modulus=1.86e9,
            stroke=70, offset=35,
        ),
        valve=RectifyingCheckValve(
            discharge_coefficient=.61, area_max=.002, area_min=1e-8,
            pressure_max=1.5e6, pressure_min=0,
            density=850, switch_gain=200,
            opening_gain=8.137375096991404e-6,
        ),
        high_accumulator=GasChargedAccumulator(8.5, 2784.7 * 6894.75),
        low_accumulator=GasChargedAccumulator(8.5, 1392.4 * 6894.75),
        motor=ConstantEfficiencyHydraulicMotor(120e-6, .9, .85),
        generator=EquivalentCircuitGenerator(.8, .8, .8, .8, .8),
        load_controller=DiscretePILoadController(1000, .001, .001),
        initial_pressure_a=2.1333e7,
        initial_pressure_b=2.1333e7,
    )


if __name__ == "__main__":
    wave = pm_equal_energy_components(
        HYDRO, significant_height=2.5, peak_period=8,
        directions=(0,), spreading=(1,), seed=1,
    )
    result = run_rm3_rectified_hydraulic(
        HYDRO, wave, build_pto(), dt=.01, end_time=20, ramp_time=5,
    )
    print(f"Samples: {len(result.time)}")
    print(f"Peak float heave: {abs(result.body_heave[:, 0]).max():.3f} m")
    print(f"Peak cylinder force: {abs(result.pto_force).max():.0f} N")

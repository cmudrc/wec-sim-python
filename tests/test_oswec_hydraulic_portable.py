"""Gate both published OSWEC hydraulic PTO layouts in ordinary Python CI."""

from hashlib import sha256
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


ROOT = Path(__file__).parent
SOURCE = ROOT / "fixtures/matlab/oswec_hydraulic/source_50ms.npz"
HYDRO = ROOT / "fixtures/matlab/oswec_multi_wave/oswec.h5"
SOURCE_SHA256 = "58c80c44eae64683da6d461b2942dd7c5fcf398fa06e48d9a49d4a389876b392"
HYDRO_SHA256 = "b9ccf3a32e0893c23b1645c15a2baef4a0c6f678b90d177d5797e3c65de4b7d3"


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


def _within(actual, expected, limit, label, failures):
    actual = np.asarray(actual)
    expected = np.asarray(expected)
    assert actual.shape == expected.shape, label
    assert np.isfinite(actual).all() and np.isfinite(expected).all(), label
    difference = np.abs(actual - expected)
    index = np.unravel_index(np.argmax(difference), difference.shape)
    error = float(difference[index])
    print(f"{label}: {error:.9g} at {index}, limit {limit:.9g}; "
          f"Python {actual[index]:.9g}, MATLAB {expected[index]:.9g}")
    if error >= limit:
        failures.append(f"{label}: {error:.9g} exceeds {limit:.9g} at {index}")


@pytest.mark.parametrize("name,fixed_crank", [
    ("adjustable", False),
    ("fixed", True),
])
def test_published_400s_hydraulic_motion_pressure_and_power(name, fixed_crank):
    failures = []

    def check(actual, expected, limit, label):
        _within(actual, expected, limit, label, failures)

    assert sha256(SOURCE.read_bytes()).hexdigest() == SOURCE_SHA256
    assert sha256(HYDRO.read_bytes()).hexdigest() == HYDRO_SHA256
    with np.load(SOURCE, allow_pickle=False) as archive:
        source = {
            key.removeprefix(f"{name}_"): archive[key]
            for key in archive.files if key.startswith(f"{name}_")
        }
    assert source["components"].shape == (500, 4)
    np.testing.assert_allclose(source["time"], np.arange(8001) * .05,
                               rtol=0, atol=1e-9)

    components = source["components"]
    sea = pm_equal_energy_components(
        HYDRO, significant_height=2.5, peak_period=8,
        directions=(0,), spreading=(1,), phase=components[:, 3:4],
    )
    for actual, expected in (
        (sea.omega, components[:, 0]),
        (sea.spectral_amplitude, components[:, 1]),
        (sea.d_omega, components[:, 2]),
    ):
        check(actual, expected, 1e-12, "published PM components")
    linkage = FixedRodCrank(3, 1.3, 5) if fixed_crank else AdjustableRodCrank(3, 1.3, 5)
    result = run_oswec_rectified_hydraulic(
        HYDRO, sea, _published_pto(fixed_crank), linkage,
    )
    assert result.time.shape == (40_001,)
    at_source = slice(None, None, 5)
    np.testing.assert_allclose(result.time[at_source], source["time"],
                               rtol=0, atol=1e-9)

    check(result.wave_elevation[at_source], source["wave"],
            1e-10, "published wave")
    check(result.excitation_force[at_source], source["excitation"],
            1e-6, "six-component excitation")
    check(result.pitch[at_source], source["pitch"],
            .002, "flap pitch")
    check(result.pitch_velocity[at_source], source["pitch_speed"],
            .001, "flap pitch speed")
    check(result.body_position[at_source], source["body_position"],
            .008, "flap center")
    check(result.body_velocity[at_source], source["body_velocity"],
            .005, "flap center speed")

    crank = source["crank"]
    check(result.crank_torque[at_source], crank[:, 0],
            75_000, "crank torque")
    stroke, jacobian = linkage.stroke_and_jacobian(crank[:, 1])
    check(result.crank_stroke[at_source], stroke,
            .008, "crank stroke")
    check(result.crank_speed[at_source], jacobian * crank[:, 2],
            .005, "crank speed")

    cylinder = source["cylinder"]
    check(result.cylinder_force[at_source], cylinder[:, 1],
            25_000, "cylinder force")
    check(result.cylinder_pressure[at_source, 0], cylinder[:, 0],
            500_000, "chamber A")
    check(result.cylinder_pressure[at_source, 1], cylinder[:, 2],
            500_000, "chamber B")
    check(result.valve_flow[at_source], source["valve"],
            .065, "valve-port flow")
    check(result.accumulator_pressure[at_source],
            source["accumulator_pressure"], 7_000,
            "high and low accumulator pressure")

    motor = source["motor"]
    check(result.shaft_speed_rpm[at_source], motor[:, 0],
            .08, "shaft speed")
    check(result.motor_torque[at_source], motor[:, 1],
            .15, "motor torque")
    check(result.motor_flow[at_source], motor[:, 2],
            2e-7, "motor flow")
    generator = source["generator"]
    check(result.generator_current[at_source], generator[:, 0],
            .2, "generator current")
    check(result.generator_voltage[at_source], generator[:, 1],
            .25, "generator voltage")

    # Integrate every Python step against MATLAB's full 10 ms work integral;
    # the archived 50 ms samples are not used to approximate the source work.
    python_work = np.trapezoid(
        result.crank_torque * result.pitch_velocity, x=result.time,
    )
    source_work = float(source["source_work"])
    assert abs(python_work - source_work) / abs(source_work) < .01
    assert not failures, "\n".join(failures)

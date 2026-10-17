"""A small RM3 sea should drive a nearly linear two-body/PTO response."""

from pathlib import Path

import numpy as np
import pytest

from wecsim import PMWave, RegularWave, WEC


HYDRO = Path(__file__).resolve().parents[1] / "examples/data/rm3.h5"
PTO_DAMPING = 1_200_000


def _run(height: float, mode: str):
    wec = WEC("RM3 small-wave scaling")
    float_body = wec.body("float", HYDRO, inertia=(0, 21_306_090.66, 0))
    spar = wec.body("spar", HYDRO, inertia=(0, 94_407_091.24, 0))
    wec.floating_joint(float_body, spar, damping=PTO_DAMPING)
    if mode == "regular":
        wave = RegularWave(height, 8)
        settings = {}
    else:
        # Keep the same seeded realization at both heights. This tests the
        # incident-sea synthesis and the convolution-radiation motion path.
        wave = PMWave(height, 8, seed=4, phase_generator="matlab",
                      frequency_count=32)
        settings = {"radiation_memory": 2}
    return wec.run(wave, dt=.1, end_time=20, ramp_time=5, **settings)


@pytest.mark.parametrize("mode", ("regular", "pm"))
def test_small_wave_amplitude_scales_rm3_motion_and_pto_work(mode):
    low = _run(.01, mode)
    high = _run(.02, mode)
    np.testing.assert_array_equal(low.time, high.time)
    np.testing.assert_allclose(high.wave_elevation, 2 * low.wave_elevation,
                               rtol=0, atol=1e-14)
    assert np.ptp(high.wave_elevation) > .001

    # The RM3 slider geometry is nonlinear, so the body response is only
    # approximately linear even though the wave elevation scales exactly.
    for name in ("float", "spar"):
        for field in ("position", "velocity"):
            small = getattr(low.bodies[name], field)
            large = getattr(high.bodies[name], field)
            if field == "position":
                small = small - small[0]
                large = large - large[0]
            assert np.max(np.abs(large - 2 * small)) < 3e-6

    small_pto = low.ptos["relative_heave"]
    large_pto = high.ptos["relative_heave"]
    for field, limit in (("stroke", 1e-6), ("velocity", 1e-6),
                         ("force", 1), ("absorbed_power", .02)):
        factor = 4 if field == "absorbed_power" else 2
        small = getattr(small_pto, field)
        large = getattr(large_pto, field)
        assert np.max(np.abs(large - factor * small)) < limit, field

    # The passive damper must absorb, never inject, power at every sample.
    for pto in (small_pto, large_pto):
        np.testing.assert_allclose(pto.absorbed_power,
                                   PTO_DAMPING * pto.velocity**2,
                                   rtol=0, atol=1e-9)
        assert np.min(pto.absorbed_power) >= 0
    low_work = np.trapezoid(small_pto.absorbed_power, x=low.time)
    high_work = np.trapezoid(large_pto.absorbed_power, x=high.time)
    assert high_work > 1
    assert abs(high_work / low_work - 4) < .001

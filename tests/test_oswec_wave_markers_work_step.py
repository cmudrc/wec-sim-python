"""Separate published OSWEC marker-case PTO work from source step sensitivity."""

from pathlib import Path

import numpy as np

from wecsim import RegularWave, WEC, WorldPoint


ROOT = Path(__file__).resolve().parents[1]
PUBLISHED = ROOT / "tests/fixtures/matlab/oswec_wave_markers.npz"
HALF_STEP = (ROOT / "tests/fixtures/matlab/oswec_wave_markers_step"
             / "half_source.npz")
HYDRO = ROOT / "examples/data/oswec.h5"


def test_full_duration_oswec_marker_work_step_separation():
    with np.load(PUBLISHED, allow_pickle=False) as published, \
            np.load(HALF_STEP, allow_pickle=False) as half:
        time = published["time"]
        np.testing.assert_array_equal(half["time"], time)
        np.testing.assert_allclose(half["wave"], published["origin"],
                                   rtol=0, atol=1e-12)
        assert float(half["simu_dt"]) == .05
        assert float(half["dt_out"]) == .1
        assert half["pto_speed"].shape == time.shape == (4001,)
        np.testing.assert_allclose(
            half["pto_force"], -12_000 * half["pto_speed"],
            rtol=0, atol=1e-5,
        )
        np.testing.assert_allclose(
            -half["pto_power"], 12_000 * half["pto_speed"]**2,
            rtol=0, atol=1e-5,
        )

        wec = WEC("OSWEC Wave Markers work audit")
        flap = wec.body("flap", HYDRO, mass=127_000,
                        inertia=(0, 1.85e6, 0))
        base = wec.body("base", HYDRO, mass=999,
                        inertia=(999, 999, 999))
        wec.fixed_hinge(flap, base, location=WorldPoint(0, 0, -10),
                        pto_location=WorldPoint(0, 0, -8.9), damping=12_000)
        result = wec.run(RegularWave(2.5, 8), dt=.1,
                         end_time=400, ramp_time=100)
        np.testing.assert_allclose(result.time, time, rtol=0, atol=1e-10)

        for label, start, stop in (
            ("ramp", 0, 100),
            ("early steady", 100, 250),
            ("late steady", 250, 400),
            ("full", 0, 400),
        ):
            interval = (time >= start) & (time <= stop)
            published_work = -np.trapezoid(
                published["pto_power"][interval], time[interval],
            )
            half_work = -np.trapezoid(
                half["pto_power"][interval], time[interval],
            )
            python_work = np.trapezoid(
                result.ptos["hinge"].absorbed_power[interval],
                time[interval],
            )
            assert min(published_work, half_work, python_work) > 0, label
            source_change = abs(half_work - published_work) / published_work
            python_gap = abs(python_work - half_work) / half_work
            # The original 1% published-work gate fails. Hold that limit for
            # the independent Python/half-step comparison; the immutable
            # MATLAB pair exposes why it cannot certify published-step work.
            assert source_change > .01, (label, source_change)
            assert python_gap < .01, (label, python_gap)

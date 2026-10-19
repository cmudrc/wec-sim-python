"""Separate published OSWEC marker-case PTO work from source step sensitivity."""

from pathlib import Path

import numpy as np

from wecsim import RegularWave, WEC, WorldPoint


ROOT = Path(__file__).resolve().parents[1]
PUBLISHED = ROOT / "tests/fixtures/matlab/oswec_wave_markers.npz"
HALF_STEP = (ROOT / "tests/fixtures/matlab/oswec_wave_markers_step"
             / "half_source.npz")
QUARTER_STEP = (ROOT / "tests/fixtures/matlab/oswec_wave_markers_step"
                / "quarter_source.npz")
EIGHTH_STEP = (ROOT / "tests/fixtures/matlab/oswec_wave_markers_step"
              / "eighth_source.npz")
# Use the BEMIO output from this application, not the different example HDF5.
HYDRO = ROOT / "tests/fixtures/matlab/oswec_wave_markers/oswec.h5"


def test_full_duration_oswec_marker_work_step_separation():
    with np.load(PUBLISHED, allow_pickle=False) as published, \
            np.load(HALF_STEP, allow_pickle=False) as half, \
            np.load(QUARTER_STEP, allow_pickle=False) as quarter, \
            np.load(EIGHTH_STEP, allow_pickle=False) as eighth:
        time = published["time"]
        for source, step in ((half, .05), (quarter, .025),
                             (eighth, .0125)):
            np.testing.assert_array_equal(source["time"], time)
            np.testing.assert_allclose(source["wave"], published["origin"],
                                       rtol=0, atol=1e-12)
            assert float(source["simu_dt"]) == step
            assert float(source["dt_out"]) == .1
            assert source["pto_speed"].shape == time.shape == (4001,)
            np.testing.assert_allclose(
                source["pto_force"], -12_000 * source["pto_speed"],
                rtol=0, atol=1e-5,
            )
            np.testing.assert_allclose(
                -source["pto_power"], 12_000 * source["pto_speed"]**2,
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
            quarter_work = -np.trapezoid(
                quarter["pto_power"][interval], time[interval],
            )
            eighth_work = -np.trapezoid(
                eighth["pto_power"][interval], time[interval],
            )
            python_work = np.trapezoid(
                result.ptos["hinge"].absorbed_power[interval],
                time[interval],
            )
            assert min(published_work, half_work, quarter_work, eighth_work,
                       python_work) > 0, label
            source_change = abs(half_work - published_work) / published_work
            half_gap = abs(python_work - half_work) / half_work
            quarter_gap = abs(python_work - quarter_work) / quarter_work
            eighth_gap = abs(python_work - eighth_work) / eighth_work
            # The successive work differences approach a factor of four on
            # each halving. The extrapolation checks the independent Python
            # result against that measured second-order source trend.
            change_ratio = ((half_work - quarter_work)
                            / (quarter_work - eighth_work))
            extrapolated_work = (4 * eighth_work - quarter_work) / 3
            # The original 1% published-work gate fails. Hold that limit for
            # each independent Python/refined-source comparison. The second
            # and third MATLAB refinements must continue toward Python.
            assert source_change > .01, (label, source_change)
            assert half_gap < .01, (label, half_gap)
            assert quarter_gap < half_gap, (label, quarter_gap, half_gap)
            assert quarter_gap < .01, (label, quarter_gap)
            assert eighth_gap < quarter_gap, (label, eighth_gap, quarter_gap)
            assert eighth_gap < .01, (label, eighth_gap)
            assert 3.5 < change_ratio < 4.5, (label, change_ratio)
            assert abs(python_work - extrapolated_work) / extrapolated_work < .001

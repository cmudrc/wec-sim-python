"""Separate published OSWEC PTO work bias from refined physical motion."""

from pathlib import Path

import numpy as np

from wecsim import PMWave, WEC, WorldPoint


ROOT = Path(__file__).resolve().parents[1]
PUBLISHED = ROOT / "tests/fixtures/matlab/oswec_irregular_wave.npz"
REFINED = ROOT / "tests/fixtures/matlab/oswec_core_step_work/source.npz"
HYDRO = ROOT / "examples/data/oswec.h5"
WINDOWS = ((0, 400), (0, 100), (100, 250), (250, 400))


def _work(time, absorbed_power, start, stop):
    window = (time >= start) & (time <= stop)
    return float(np.trapezoid(absorbed_power[window], time[window]))


def test_matlab_oswec_fixed_step_changes_published_pto_work():
    with np.load(PUBLISHED, allow_pickle=False) as original, \
            np.load(REFINED, allow_pickle=False) as saved:
        half_pto = saved["half_pto"]
        quarter_pto = saved["quarter_pto"]
        half_wave = saved["half_wave"]
        quarter_wave = saved["quarter_wave"]
        assert half_pto.shape == (8001, 4)
        assert quarter_pto.shape == (16001, 4)
        np.testing.assert_allclose(half_pto[:, 0], np.arange(8001) * .05,
                                   rtol=0, atol=3e-11)
        np.testing.assert_allclose(quarter_pto[:, 0], np.arange(16001) * .025,
                                   rtol=0, atol=3e-11)
        np.testing.assert_array_equal(half_wave[::2, 1],
                                      original["wave_elevation"])
        np.testing.assert_array_equal(quarter_wave[::4, 1],
                                      original["wave_elevation"])
        np.testing.assert_array_equal(half_wave[::2], quarter_wave[::4])
        for trace in (half_pto, quarter_pto):
            # MATLAB reports mechanical power with the opposite sign from
            # absorbed power. The published linear PTO law is unchanged.
            np.testing.assert_allclose(trace[:, 2], -12_000 * trace[:, 1],
                                       rtol=0, atol=1e-8)
            np.testing.assert_allclose(trace[:, 3],
                                       -12_000 * trace[:, 1] ** 2,
                                       rtol=0, atol=1e-8)

        for start, stop in WINDOWS:
            coarse_work = _work(original["time"], -original["pto_power"],
                                start, stop)
            half_work = _work(half_pto[:, 0], -half_pto[:, 3], start, stop)
            quarter_work = _work(quarter_pto[:, 0], -quarter_pto[:, 3],
                                 start, stop)
            assert coarse_work > 1_000
            first_change = (coarse_work - half_work) / coarse_work
            second_change = (half_work - quarter_work) / half_work
            assert .005 < first_change < .01
            assert .001 < second_change < .003
            assert second_change < first_change / 2


def test_independent_oswec_motion_and_pto_work_against_refined_source(tmp_path):
    with np.load(PUBLISHED, allow_pickle=False) as original, \
            np.load(REFINED, allow_pickle=False) as saved:
        phases = original["wave_components"][:, 3:].copy()
        directions = original["wave_directions"].copy()
        source_wave = saved["quarter_wave"].copy()
        source_body = saved["quarter_body_active"].copy()
        source_pto = saved["quarter_pto"].copy()
    phase_file = tmp_path / "phase.csv"
    np.savetxt(phase_file, phases, delimiter=",")
    wec = WEC("Published OSWEC step-refinement audit")
    flap = wec.body("flap", HYDRO, mass=127_000, inertia=(1.85e6,) * 3)
    base = wec.body("base", HYDRO, mass=999, inertia=(999,) * 3)
    wec.fixed_hinge(flap, base, location=WorldPoint(0, 0, -10),
                    pto_location=WorldPoint(0, 0, -8.9), damping=12_000)
    result = wec.run(
        PMWave(2.5, 8, directions=tuple(directions[:, 0]),
               spreading=tuple(directions[:, 1]), phase_file=phase_file),
        dt=.025, end_time=400, ramp_time=100, radiation_memory=30,
    )
    np.testing.assert_allclose(result.time, source_body[:, 0], rtol=0,
                               atol=3e-11)
    np.testing.assert_allclose(result.wave_elevation, source_wave[:, 1],
                               rtol=0, atol=1e-11)
    motion = result.bodies["flap"]
    for actual, expected, limits in (
        (motion.position[:, [0, 2, 4]], source_body[:, 1:4],
         (.002, .0005, .0005)),
        (motion.velocity[:, [0, 2, 4]], source_body[:, 4:7],
         (.003, .001, .001)),
    ):
        difference = np.max(np.abs(actual - expected), axis=0)
        assert np.all(difference < limits), (difference, limits)
    pto = result.ptos["hinge"]
    assert np.max(np.abs(pto.velocity - source_pto[:, 1])) < .001
    assert np.max(np.abs(pto.force - source_pto[:, 2])) < 5
    for start, stop in WINDOWS:
        source_work = _work(source_pto[:, 0], -source_pto[:, 3], start, stop)
        python_work = _work(result.time, pto.absorbed_power, start, stop)
        assert source_work > 1_000
        relative_gap = abs(python_work - source_work) / source_work
        assert relative_gap < .001, (
            f"OSWEC PTO work {start}–{stop} s: {relative_gap:.4%}"
        )

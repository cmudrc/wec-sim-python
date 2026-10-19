"""Check converged core OSWEC energy without promoting its coarse source step."""

from pathlib import Path

import numpy as np

from wecsim import PMWave, WEC, WorldPoint


ROOT = Path(__file__).resolve().parents[1]
PUBLISHED = ROOT / "tests/fixtures/matlab/oswec_irregular_wave.npz"
REFINED = ROOT / "tests/fixtures/matlab/oswec_core_step_work/source.npz"
EIGHTH = ROOT / "tests/fixtures/matlab/oswec_core_step_work/eighth_source.npz"
HYDRO = ROOT / "examples/data/oswec.h5"
WINDOWS = (("full", 0, 400), ("ramp", 0, 100),
           ("early steady", 100, 250), ("late steady", 250, 400))


def _work(time, power, start, stop):
    selected = (time >= start) & (time <= stop)
    return float(np.trapezoid(power[selected], time[selected]))


def test_core_oswec_eighth_step_motion_and_extrapolated_work(tmp_path):
    with np.load(PUBLISHED, allow_pickle=False) as published, \
            np.load(REFINED, allow_pickle=False) as refined, \
            np.load(EIGHTH, allow_pickle=False) as eighth:
        phases = published["wave_components"][:, 3:].copy()
        directions = published["wave_directions"].copy()
        half_pto = refined["half_pto"].copy()
        quarter_pto = refined["quarter_pto"].copy()
        source_pto = eighth["pto"].copy()
        source_wave = eighth["wave"].copy()
        source_body = eighth["body_active"].copy()
        assert float(eighth["simu_dt"]) == .0125
        assert source_pto.shape == (32001, 4)
        assert source_wave.shape == (32001, 2)
        assert source_body.shape == (32001, 7)
        np.testing.assert_allclose(source_pto[:, 0],
                                   np.arange(32001) * .0125,
                                   rtol=0, atol=3e-11)
        np.testing.assert_array_equal(source_wave[::8, 1],
                                      published["wave_elevation"])
        np.testing.assert_array_equal(source_wave[:, 0], source_pto[:, 0])
        np.testing.assert_array_equal(source_body[:, 0], source_pto[:, 0])
        np.testing.assert_allclose(source_pto[:, 2],
                                   -12_000 * source_pto[:, 1],
                                   rtol=0, atol=1e-8)
        np.testing.assert_allclose(source_pto[:, 3],
                                   -12_000 * source_pto[:, 1] ** 2,
                                   rtol=0, atol=1e-8)

    phase_file = tmp_path / "phase.csv"
    np.savetxt(phase_file, phases, delimiter=",")
    python_work = {}
    for dt in (.025, .0125, .00625):
        wec = WEC("Core OSWEC convergence audit")
        flap = wec.body("flap", HYDRO, mass=127_000,
                        inertia=(1.85e6,) * 3)
        base = wec.body("base", HYDRO, mass=999, inertia=(999,) * 3)
        wec.fixed_hinge(flap, base, location=WorldPoint(0, 0, -10),
                        pto_location=WorldPoint(0, 0, -8.9), damping=12_000)
        result = wec.run(
            PMWave(2.5, 8, directions=tuple(directions[:, 0]),
                   spreading=tuple(directions[:, 1]), phase_file=phase_file),
            dt=dt, end_time=400, ramp_time=100, radiation_memory=30,
        )
        pto = result.ptos["hinge"]
        python_work[dt] = {
            label: _work(result.time, pto.absorbed_power, start, stop)
            for label, start, stop in WINDOWS
        }
        if dt == .0125:
            np.testing.assert_allclose(result.time, source_pto[:, 0],
                                       rtol=0, atol=3e-11)
            np.testing.assert_allclose(result.wave_elevation,
                                       source_wave[:, 1], rtol=0, atol=1e-11)
            motion = result.bodies["flap"]
            for actual, expected, limits in (
                (motion.position[:, [0, 2, 4]], source_body[:, 1:4],
                 (.002, .0005, .0005)),
                (motion.velocity[:, [0, 2, 4]], source_body[:, 4:7],
                 (.003, .001, .001)),
            ):
                difference = np.max(np.abs(actual - expected), axis=0)
                assert np.all(difference < limits), (difference, limits)
            assert np.max(np.abs(pto.velocity - source_pto[:, 1])) < .001
            assert np.max(np.abs(pto.force - source_pto[:, 2])) < 5

    for label, start, stop in WINDOWS:
        half = _work(half_pto[:, 0], -half_pto[:, 3], start, stop)
        quarter = _work(quarter_pto[:, 0], -quarter_pto[:, 3], start, stop)
        eighth = _work(source_pto[:, 0], -source_pto[:, 3], start, stop)
        python_quarter = python_work[.025][label]
        python_eighth = python_work[.0125][label]
        python_sixteenth = python_work[.00625][label]
        assert min(half, quarter, eighth, python_quarter, python_eighth,
                   python_sixteenth) > 0, label
        # Keep the earlier 0.1% derived-step work limit at the finer step.
        assert abs(python_eighth - eighth) / eighth < .001, label
        source_ratio = (half - quarter) / (quarter - eighth)
        python_ratio = ((python_eighth - python_quarter)
                        / (python_sixteenth - python_eighth))
        assert 3.5 < source_ratio < 4.5, (label, source_ratio)
        assert 3.5 < python_ratio < 4.5, (label, python_ratio)
        source_extrapolated = (4 * eighth - quarter) / 3
        python_extrapolated = (4 * python_sixteenth - python_eighth) / 3
        # This derived convergence regression does not pass the published
        # 0.1 s accumulated-work comparison or change the case badge.
        assert (abs(python_extrapolated - source_extrapolated)
                / source_extrapolated < 1e-4), label

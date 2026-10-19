"""Separate OSWEC's ode4 step from its published wave and force clocks."""

from pathlib import Path

import numpy as np

from wecsim import PMWave, WEC, WorldPoint


ROOT = Path(__file__).resolve().parents[1]
PUBLISHED = ROOT / "tests/fixtures/matlab/oswec_irregular_wave.npz"
SOLVER_ONLY = ROOT / "tests/fixtures/matlab/oswec_core_fixedstep_only/source.npz"
HYDRO = ROOT / "examples/data/oswec.h5"
WINDOWS = ((0, 400), (0, 100), (100, 250), (250, 400))


def _work(time, absorbed_power, start, stop):
    selected = (time >= start - 1e-9) & (time <= stop + 1e-9)
    return float(np.trapezoid(absorbed_power[selected], time[selected]))


def test_matlab_solver_step_alone_changes_published_oswec_work():
    with np.load(PUBLISHED, allow_pickle=False) as original, \
            np.load(SOLVER_ONLY, allow_pickle=False) as source:
        for label, step in (("coarse", .1), ("half", .05), ("quarter", .025)):
            pto, wave = source[f"{label}_pto"], source[f"{label}_wave"]
            settings = source[f"{label}_settings"]
            assert pto.shape == (4001, 4)
            assert wave.shape == (4001, 2)
            assert source[f"{label}_pitch"].shape == (4001,)
            assert np.isfinite(pto).all()
            np.testing.assert_allclose(settings, [.1, .1, .1, .1, step],
                                       rtol=0, atol=1e-12)
            np.testing.assert_allclose(pto[:, 0], original["time"],
                                       rtol=0, atol=1e-10)
            np.testing.assert_array_equal(wave[:, 1],
                                          original["wave_elevation"])
            np.testing.assert_allclose(pto[:, 2], -12_000 * pto[:, 1],
                                       rtol=0, atol=1e-8)
            np.testing.assert_allclose(pto[:, 3], -12_000 * pto[:, 1] ** 2,
                                       rtol=0, atol=1e-8)

        # The fresh local 0.1 s run reproduces the pre-existing source record.
        np.testing.assert_allclose(source["coarse_pto"][:, 1],
                                   original["pto_speed"], rtol=0, atol=2e-14)
        np.testing.assert_allclose(source["coarse_pto"][:, 3],
                                   original["pto_power"], rtol=0, atol=1e-8)
        np.testing.assert_allclose(source["coarse_pitch"],
                                   original["body_position"][0, :, 4],
                                   rtol=0, atol=2e-13)

        for start, stop in WINDOWS:
            work = [
                _work(source[f"{label}_pto"][:, 0],
                      -source[f"{label}_pto"][:, 3], start, stop)
                for label in ("coarse", "half", "quarter")
            ]
            assert min(work) > 1_000
            first_change = (work[0] - work[1]) / work[0]
            second_change = (work[1] - work[2]) / work[1]
            # All source force/output clocks remain at 0.1 s. This is a
            # source self-change diagnostic, not a published-case parity gate.
            assert .0085 < first_change < .0115
            assert .0018 < second_change < .0025
            assert second_change < first_change / 3


def test_independent_oswec_work_against_solver_refined_source(tmp_path):
    with np.load(PUBLISHED, allow_pickle=False) as original, \
            np.load(SOLVER_ONLY, allow_pickle=False) as source:
        phases = original["wave_components"][:, 3:].copy()
        directions = original["wave_directions"].copy()
        quarter_wave = source["quarter_wave"][:, 1].copy()
        quarter_active = source["quarter_active"].copy()
        quarter_pto = source["quarter_pto"].copy()
    phase_file = tmp_path / "phase.csv"
    np.savetxt(phase_file, phases, delimiter=",")
    wec = WEC("OSWEC solver-step isolation")
    flap = wec.body("flap", HYDRO, mass=127_000, inertia=(1.85e6,) * 3)
    base = wec.body("base", HYDRO, mass=999, inertia=(999,) * 3)
    wec.fixed_hinge(flap, base, location=WorldPoint(0, 0, -10),
                    pto_location=WorldPoint(0, 0, -8.9), damping=12_000)
    result = wec.run(
        PMWave(2.5, 8, directions=tuple(directions[:, 0]),
               spreading=tuple(directions[:, 1]), phase_file=phase_file),
        dt=.1, end_time=400, ramp_time=100, radiation_memory=30,
    )
    np.testing.assert_allclose(result.wave_elevation, quarter_wave,
                               rtol=0, atol=1e-11)
    motion = result.bodies["flap"]
    for actual, expected, limits in (
        (motion.position[:, [0, 2, 4]], quarter_active[:, :3],
         (.03, .005, .005)),
        (motion.velocity[:, [0, 2, 4]], quarter_active[:, 3:],
         (.015, .004, .003)),
    ):
        differences = np.max(np.abs(actual - expected), axis=0)
        assert np.all(differences < limits), (differences, limits)
    pto = result.ptos["hinge"]
    assert np.max(np.abs(pto.velocity - quarter_pto[:, 1])) < .003
    assert np.max(np.abs(pto.force - quarter_pto[:, 2])) < 30
    assert np.max(np.abs(pto.absorbed_power + quarter_pto[:, 3])) < 7
    for start, stop in WINDOWS:
        source_work = _work(quarter_pto[:, 0], -quarter_pto[:, 3],
                            start, stop)
        python_work = _work(result.time, pto.absorbed_power, start, stop)
        assert abs(python_work - source_work) / source_work < .002

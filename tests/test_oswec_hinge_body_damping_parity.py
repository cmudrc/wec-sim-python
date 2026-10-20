"""Pair a configured hinged flap body damper with a local MATLAB source run."""

from pathlib import Path

import numpy as np

from wecsim import RegularWave, WEC, WorldPoint


ROOT = Path(__file__).resolve().parents[1]
HYDRO = ROOT / "tests/fixtures/matlab/oswec_multi_wave/oswec.h5"
SOURCE = ROOT / "tests/fixtures/matlab/oswec_hinge_body_damping/source.npz"
PUBLISHED = ROOT / "tests/fixtures/matlab/oswec_nonhydro/source.npz"
AXES = (0, 2, 4)
POSITION_LIMITS = (0.03, 0.012, 0.007)
SPEED_LIMITS = (0.025, 0.01, 0.006)
DELTA_POSITION_LIMITS = (0.001, 0.0005, 0.00015)
DELTA_SPEED_LIMITS = (0.001, 0.0005, 0.00025)


def _maximum(left, right):
    assert left.shape == right.shape
    return float(np.max(np.abs(left - right)))


def _run(matrix, *, source_ode4=False):
    wec = WEC("OSWEC hinged flap with body damping")
    flap = wec.body("flap", HYDRO, mass=127_000,
                    inertia=(0, 1.85e6, 0), linear_damping=matrix)
    base = wec.fixed_body("base", center_gravity=(0, 0, -10.9),
                          mass=999, inertia=(1, 1, 1))
    wec.fixed_hinge(flap, base, location=WorldPoint(0, 0, -10),
                    pto_location=WorldPoint(0, 0, -8.9),
                    added_mass_scheme=("simulink_delay" if source_ode4
                                       else "implicit"))
    return wec.run(RegularWave(2.5, 8), dt=.1, end_time=150,
                   ramp_time=100,
                   **({"solver": "wecsim_ode4"} if source_ode4 else {}))


def test_body_damping_matches_source_and_changes_hinge_motion():
    with np.load(SOURCE, allow_pickle=False) as saved:
        flap = saved["flap"].copy()
        base = saved["base"].copy()
        pto = saved["pto"].copy()
        fine_flap = saved["fine_flap"].copy()
        fine_base = saved["fine_base"].copy()
        fine_pto = saved["fine_pto"].copy()
        source_linear_force = saved["linear_force"].copy()
        matrix = saved["damping"].copy()
        fixed_steps = saved["fixed_steps"].copy()
    assert flap.shape == base.shape == (2, 1501, 13)
    assert pto.shape == (2, 1501, 19)
    assert fine_flap.shape == fine_base.shape == flap.shape
    assert fine_pto.shape == pto.shape
    assert source_linear_force.shape == (2, 1501, 6)
    np.testing.assert_array_equal(fixed_steps, [.025, .025])
    np.testing.assert_array_equal(matrix, np.diag([100_000, 0, 0, 0, 2_000_000, 0]))
    with np.load(PUBLISHED, allow_pickle=False) as published:
        # The unchanged arm verifies that the local source copy reproduces
        # the pinned published case before testing the derived setting.
        assert _maximum(flap[0], published["flap"][:1501, :13]) < 1e-7
        assert _maximum(base[0], published["base"][:1501, :13]) < 1e-7
    for timeline in flap[:, :, 0]:
        np.testing.assert_allclose(timeline, np.arange(1501) * .1,
                                   rtol=0, atol=1e-10)
    np.testing.assert_allclose(base[0], base[1], rtol=0, atol=1e-10)
    # A body damping coefficient, separate from the zero-damping PTO, must
    # visibly affect the source trajectory before a paired gate is useful.
    assert _maximum(flap[0, :, 5], flap[1, :, 5]) > .02
    assert np.max(np.abs(pto[:, :, 17])) < 1e-6
    # The source logs the damping term as +D v; its force balance subtracts
    # that term. This checks all six components on source-generated motion.
    assert np.max(np.abs(source_linear_force[0])) == 0
    assert _maximum(source_linear_force[1],
                    flap[1, :, 7:13] @ matrix.T) < 1e-8

    control, damped = _run(None), _run(matrix)
    for index, result in enumerate((control, damped)):
        np.testing.assert_allclose(result.time, flap[index, :, 0],
                                   rtol=0, atol=1e-10)
        assert result.case["bodies"][0].get("linear_damping") == (
            None if index == 0 else matrix.tolist()
        )
        actual = result.bodies["flap"]
        for axis, position_limit, speed_limit in zip(
                AXES, POSITION_LIMITS, SPEED_LIMITS, strict=True):
            assert _maximum(actual.position[:, axis],
                            flap[index, :, 1 + axis]) < position_limit
            assert _maximum(actual.velocity[:, axis],
                            flap[index, :, 7 + axis]) < speed_limit
        assert _maximum(result.bodies["base"].position,
                        base[index, :, 1:7]) < 1e-10
        assert _maximum(result.bodies["base"].velocity,
                        base[index, :, 7:13]) < 1e-10
        assert _maximum(result.ptos["hinge"].force,
                        pto[index, :, 13 + 4]) < 1e-6

    # Solver-only refinement retains the 0.1 s WEC-Sim wave and output
    # clocks. The unchanged differential limits failed against the coarse
    # 0.1 s MATLAB arm; retain that as a measured limitation below.
    assert _maximum(base, fine_base) < 1e-10
    assert np.max(np.abs(fine_pto[:, :, 17])) < 1e-6
    for axis, position_limit, speed_limit in zip(
            AXES, DELTA_POSITION_LIMITS, DELTA_SPEED_LIMITS, strict=True):
        source_position_change = (fine_flap[1, :, 1 + axis]
                                  - fine_flap[0, :, 1 + axis])
        python_position_change = (damped.bodies["flap"].position[:, axis]
                                  - control.bodies["flap"].position[:, axis])
        assert _maximum(python_position_change,
                        source_position_change) < position_limit
        source_speed_change = (fine_flap[1, :, 7 + axis]
                               - fine_flap[0, :, 7 + axis])
        python_speed_change = (damped.bodies["flap"].velocity[:, axis]
                               - control.bodies["flap"].velocity[:, axis])
        assert _maximum(python_speed_change, source_speed_change) < speed_limit
    coarse_surge_change = flap[1, :, 1] - flap[0, :, 1]
    python_surge_change = (damped.bodies["flap"].position[:, 0]
                           - control.bodies["flap"].position[:, 0])
    assert _maximum(python_surge_change, coarse_surge_change) > DELTA_POSITION_LIMITS[0]


def test_opt_in_source_step_pairs_both_damping_variants():
    with np.load(SOURCE, allow_pickle=False) as saved:
        flap = saved["flap"].copy()
        matrix = saved["damping"].copy()
    for index, damping in enumerate((None, matrix)):
        result = _run(damping, source_ode4=True)
        actual = result.bodies["flap"]
        for axis in AXES:
            assert _maximum(actual.position[:, axis],
                            flap[index, :, 1 + axis]) < 1e-9
            assert _maximum(actual.velocity[:, axis],
                            flap[index, :, 7 + axis]) < 1e-9

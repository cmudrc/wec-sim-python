"""Gate the opt-in published-step OSWEC nonlinear hinge motion."""

from pathlib import Path

import numpy as np
import pytest

from wecsim import RegularWave, WEC, WorldPoint


ROOT = Path(__file__).parent
SOURCE = ROOT / "fixtures/matlab/oswec_nonlinear_viz_published_motion.npz"
# These existing fixtures have the pinned regular-frequency coefficients and
# byte-identical flap STL; no source pose or force trace is read by the solve.
HYDRO = ROOT / "fixtures/matlab/oswec_multi_wave/oswec.h5"
FLAP = ROOT / "test_objects/test_bodyclass/testData/body_3_test/flap.stl"


def _wec(*, delayed):
    wec = WEC("Published nonlinear OSWEC")
    flap = wec.body(
        "flap", HYDRO, mass=127000, inertia=(1.85e6,) * 3,
        geometry_file=FLAP, nonlinear_hydro="instantaneous",
    )
    base = wec.fixed_body(
        "base", center_gravity=(0, 0, -10.9),
        mass=999, inertia=(1, 1, 1),
    )
    wec.fixed_hinge(
        flap, base, location=WorldPoint(0, 0, -10),
        pto_location=WorldPoint(0, 0, -8.9),
        added_mass_scheme="simulink_delay" if delayed else "implicit",
    )
    return wec


def test_published_nonlinear_motion_through_public_source_solver():
    with np.load(SOURCE, allow_pickle=False) as source:
        time = source["time"]
        pose = source["pose"]
        speed = source["velocity"]
    assert time.shape == (1201,)
    assert pose.shape == speed.shape == (1201, 6)

    wave = RegularWave(2.5, 8)
    result = _wec(delayed=True).run(
        wave, dt=.1, end_time=120, ramp_time=40,
        solver="wecsim_ode4",
    )
    np.testing.assert_allclose(result.time, time, rtol=0, atol=1e-12)
    position = result.bodies["flap"].position
    velocity = result.bodies["flap"].velocity
    for axis, position_limit, velocity_limit in (
        (0, .002, .002),     # surge: m, m/s
        (2, .0008, .001),    # heave: m, m/s
        (4, .00045, .0005),  # pitch: rad, rad/s
    ):
        assert np.max(np.abs(position[:, axis] - pose[:, axis])) < position_limit
        assert np.max(np.abs(velocity[:, axis] - speed[:, axis])) < velocity_limit
    np.testing.assert_array_equal(result.ptos["hinge"].force,
                                  np.zeros_like(result.time))

    # The physical default remains distinct from this step-sensitive source
    # trajectory. The source solver is an opt-in numerical reproduction.
    physical = _wec(delayed=False).run(
        wave, dt=.1, end_time=120, ramp_time=40,
    )
    assert np.max(np.abs(physical.bodies["flap"].position[:, 0]
                         - pose[:, 0])) > .03
    assert np.max(np.abs(physical.bodies["flap"].position[:, 4]
                         - pose[:, 4])) > .006


def test_regular_nonlinear_source_solver_needs_both_opt_ins():
    wave = RegularWave(2.5, 8)
    with pytest.raises(ValueError, match="simulink_delay"):
        _wec(delayed=True).to_case(wave, dt=.1, end_time=1)
    with pytest.raises(ValueError, match="source ode4"):
        _wec(delayed=False).to_case(
            wave, dt=.1, end_time=1, solver="wecsim_ode4",
        )

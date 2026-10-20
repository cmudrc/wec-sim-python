"""Separate published RM3 End Stops force-law parity from solver scheduling.

The second gate prescribes only MATLAB's accepted major-step times. It checks
Python motion and work, but cannot qualify independent published-case parity.
"""

from pathlib import Path

import numpy as np

from wecsim import LinearHardStops
from wecsim.dynamics.rm3Regular import solve_rm3_regular


ROOT = Path(__file__).resolve().parents[1]
HYDRO = ROOT / "examples/data/rm3_end_stops.h5"
FIXTURES = ROOT / "tests/fixtures/matlab"
AXES = (0, 2, 4)


def _solve(end_time, *, major_times=None):
    return solve_rm3_regular(
        HYDRO, pto_hard_stops=LinearHardStops(-.6, .6, 1e8, 1e8),
        added_mass_scheme="simulink_delay", solver="wecsim_ode45",
        _source_major_times=major_times,
        dt=.1, end_time=end_time, ramp_time=100,
    )


def _maximum(actual, expected, limit, label):
    assert actual.shape == expected.shape, label
    error = float(np.max(np.abs(actual - expected)))
    assert np.isfinite(error) and error < limit, (
        f"{label}: {error:.9g} >= {limit:.9g}"
    )


def test_independently_selected_published_prefix_through_first_contact():
    """No source states, forces, or solver times enter the 84.3 s solve."""
    result = _solve(84.3)
    with np.load(FIXTURES / "rm3_end_stops_published_prefix.npz") as source:
        length = len(result.time)
        np.testing.assert_allclose(result.time, source["time"][:length],
                                   rtol=0, atol=1e-10)
        for kind, actual, limits in (
                ("position", result.body_position[:, :, AXES],
                 (1e-5, 1e-5, 1e-6)),
                ("velocity", result.body_velocity[:, :, AXES],
                 (1e-5, 1e-5, 1e-6))):
            expected = source[f"body_{kind}"][:length]
            for body in range(2):
                for axis, limit in enumerate(limits):
                    _maximum(actual[:, body, axis], expected[:, body, axis],
                             limit, f"{kind} body {body + 1} axis {AXES[axis]}")
        for label, actual, axis, limit in (
                ("stroke", result.pto_stroke, 0, 1e-6),
                ("speed", result.pto_velocity, 1, 1e-5),
                ("force", result.pto_force, 2, 10)):
            _maximum(actual, source["pto"][:length, axis], limit, label)


def test_published_full_duration_with_source_major_step_schedule():
    """Gate the force law while explicitly retaining the timing condition."""
    with np.load(FIXTURES / "rm3_end_stops_internal_steps.npz") as steps:
        major_times = steps["internal_times"]
    result = _solve(400, major_times=major_times)
    with np.load(FIXTURES / "rm3_end_stops_published_full.npz") as source, \
         np.load(FIXTURES / "rm3_end_stops_maxstep_400s.npz") as archived, \
         np.load(FIXTURES / "rm3_end_stops_published_prefix.npz") as prefix:
        np.testing.assert_allclose(source["time"], result.time,
                                   rtol=0, atol=1e-10)
        np.testing.assert_array_equal(source["wave"], archived["wave"][:, 1])
        np.testing.assert_array_equal(source["pto"], archived["published_pto"])
        for channel in ("body_position", "body_velocity", "pto"):
            np.testing.assert_array_equal(source[channel][:1201], prefix[channel])
        for kind, actual, limits in (
                ("position", result.body_position[:, :, AXES],
                 (.01, .01, 1e-4)),
                ("velocity", result.body_velocity[:, :, AXES],
                 (.05, .05, .005))):
            expected = source[f"body_{kind}"]
            for body in range(2):
                for axis, limit in enumerate(limits):
                    _maximum(actual[:, body, axis], expected[:, body, axis],
                             limit, f"{kind} body {body + 1} axis {AXES[axis]}")
        for label, actual, axis, limit in (
                ("stroke", result.pto_stroke, 0, .01),
                ("speed", result.pto_velocity, 1, .05),
                ("force", result.pto_force, 2, 500_000)):
            _maximum(actual, source["pto"][:, axis], limit, label)
        python_work = np.trapezoid(1_200_000 * result.pto_velocity**2,
                                   result.time)
        matlab_work = np.trapezoid(1_200_000 * source["pto"][:, 1]**2,
                                   source["time"])
        assert abs(python_work - matlab_work) / matlab_work < .001

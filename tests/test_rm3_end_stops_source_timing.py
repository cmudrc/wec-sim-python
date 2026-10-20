"""Separate published RM3 End Stops force-law parity from solver scheduling.

The second gate prescribes only MATLAB's accepted major-step times. It checks
Python motion and work, but cannot qualify independent published-case parity.
"""

from pathlib import Path

import numpy as np

from wecsim import LinearHardStops, RegularWave, WEC
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


def test_independently_selected_source_state_steps_through_120_seconds():
    """Advance through source off-grid steps without supplying their times."""
    result = _solve(120)
    with np.load(FIXTURES / "rm3_end_stops_published_full.npz") as source:
        length = len(result.time)
        np.testing.assert_allclose(result.time, source["time"][:length],
                                   rtol=0, atol=1e-10)
        for kind, actual, limits in (
                ("position", result.body_position[:, :, AXES],
                 (5e-5, 1e-4, 1e-6)),
                ("velocity", result.body_velocity[:, :, AXES],
                 (1e-4, 5e-4, 1e-6))):
            expected = source[f"body_{kind}"][:length]
            for body in range(2):
                for axis, limit in enumerate(limits):
                    _maximum(actual[:, body, axis], expected[:, body, axis],
                             limit, f"{kind} body {body + 1} axis {AXES[axis]}")
        for label, actual, axis, limit in (
                ("stroke", result.pto_stroke, 0, 1e-4),
                ("speed", result.pto_velocity, 1, 5e-4),
                ("force", result.pto_force, 2, 5_000)):
            _maximum(actual, source["pto"][:length, axis], limit, label)
        # The first off-grid steps follow 84.3 s. This narrower 90 s gate
        # catches a regression to fixed .1 s steps before the 120 s limit.
        for label, actual, axis, limit in (
                ("stroke", result.pto_stroke, 0, 1e-6),
                ("speed", result.pto_velocity, 1, 1e-5),
                ("force", result.pto_force, 2, 10)):
            _maximum(actual[:901], source["pto"][:901, axis], limit,
                     f"90 s {label}")
        python_work = np.trapezoid(1_200_000 * result.pto_velocity**2,
                                   result.time)
        matlab_work = np.trapezoid(1_200_000 * source["pto"][:length, 1]**2,
                                   source["time"][:length])
        assert abs(python_work - matlab_work) / matlab_work < 1e-5


def _public_end_stops_wec():
    """Build the published RM3 configuration without source motion inputs."""
    wec = WEC("RM3 End Stops")
    float_body = wec.body("float", HYDRO, inertia=(0, 21_306_090.66, 0))
    spar_body = wec.body("spar", HYDRO, inertia=(0, 94_407_091.24, 0))
    wec.floating_joint(
        float_body, spar_body, damping=1_200_000,
        hard_stops=LinearHardStops(-.6, .6, 1e8, 1e8),
        added_mass_scheme="simulink_delay",
    )
    return wec


def test_public_wec_run_pairs_published_end_stops_through_120_seconds():
    """The configured public WEC uses its own adaptive steps and wave input."""
    result = _public_end_stops_wec().run(
        RegularWave(2.5, 8), dt=.1, end_time=120,
        ramp_time=100, solver="wecsim_ode45",
    )

    with np.load(FIXTURES / "rm3_end_stops_published_full.npz") as source:
        length = len(result.time)
        np.testing.assert_allclose(result.time, source["time"][:length],
                                   rtol=0, atol=1e-10)
        for kind, limits in (
                ("position", (5e-5, 1e-4, 1e-6)),
                ("velocity", (1e-4, 5e-4, 1e-6))):
            actual = np.stack([
                getattr(result.bodies[name], kind)[:, AXES]
                for name in ("float", "spar")
            ], axis=1)
            expected = source[f"body_{kind}"][:length]
            for body in range(2):
                for axis, limit in enumerate(limits):
                    _maximum(actual[:, body, axis], expected[:, body, axis],
                             limit, f"public {kind} body {body + 1} axis {AXES[axis]}")
        pto = result.ptos["relative_heave"]
        for label, actual, axis, limit in (
                ("stroke", pto.stroke, 0, 1e-4),
                ("speed", pto.velocity, 1, 5e-4),
                ("force", pto.force, 2, 5_000)):
            _maximum(actual, source["pto"][:length, axis], limit,
                     f"public {label}")
        python_work = np.trapezoid(1_200_000 * pto.velocity**2, result.time)
        matlab_work = np.trapezoid(
            1_200_000 * source["pto"][:length, 1]**2,
            source["time"][:length],
        )
        assert abs(python_work - matlab_work) / matlab_work < 1e-5


def test_public_wec_run_published_400s_damper_work():
    """Gate independent full-run PTO work while motion remains qualified."""
    result = _public_end_stops_wec().run(
        RegularWave(2.5, 8), dt=.1, end_time=400,
        ramp_time=100, solver="wecsim_ode45",
    )
    with np.load(FIXTURES / "rm3_end_stops_published_full.npz") as source:
        np.testing.assert_allclose(result.time, source["time"],
                                   rtol=0, atol=1e-10)
        source_work = np.trapezoid(
            1_200_000 * source["pto"][:, 1]**2, source["time"],
        )
        pto = result.ptos["relative_heave"]
        np.testing.assert_allclose(
            pto.absorbed_power, 1_200_000 * pto.velocity**2,
            rtol=0, atol=1e-8,
        )
        python_work = np.trapezoid(pto.absorbed_power, result.time)
        assert source_work > 0
        # Reuse the existing 0.7% full-run End Stops damper-work limit.
        gap = abs(python_work - source_work) / source_work
        assert gap < .007, f"published 400 s damper work gap: {gap:.4%}"


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

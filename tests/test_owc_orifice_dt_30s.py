"""Keep the OWC source clock audit separate from published-case parity."""

from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid


FIXTURES = Path(__file__).parent / "fixtures/matlab"


def test_owc_source_surge_does_not_converge_monotonically_over_30s():
    with np.load(FIXTURES / "owc_orifice_dt_30s/source.npz",
                 allow_pickle=False) as saved, \
            np.load(FIXTURES / "owc_orifice_full_maxstep/source.npz",
                    allow_pickle=False) as full:
        time = saved["time"]
        interval = saved["interval"]
        motion = saved["motion"]
        velocity = saved["velocity"]
        wave = saved["wave"]
        components = saved["components"]
        work_internal = saved["work_internal"]
        work_output = saved["work_output"]

        assert time.shape == (6001,)
        np.testing.assert_allclose(time, np.arange(6001) * .005,
                                   rtol=0, atol=1e-11)
        np.testing.assert_array_equal(interval, [.005, .0025, .00125])
        assert motion.shape == velocity.shape == (3, 6001, 4)
        assert wave.shape == (3, 6001)
        assert components.shape == (3, 500, 4)
        assert np.isfinite(motion).all() and np.isfinite(velocity).all()

        # The first arm is the previously archived 130 s published-clock
        # source, not a new short-case baseline with a different startup.
        np.testing.assert_array_equal(motion[0], full["nominal_motion"][:6001])
        np.testing.assert_array_equal(wave[0], full["nominal_wave"][:6001])
        np.testing.assert_array_equal(components[0],
                                      full["nominal_components"])
        np.testing.assert_array_equal(
            wave[1:], np.broadcast_to(wave[0], wave[1:].shape))
        np.testing.assert_array_equal(
            components[1:],
            np.broadcast_to(components[0], components[1:].shape))

    first_surge_change = float(np.max(np.abs(motion[1, :, 0]
                                             - motion[0, :, 0])))
    second_surge_change = float(np.max(np.abs(motion[2, :, 0]
                                              - motion[1, :, 0])))
    # The smaller source step makes the surge self-change larger. This is a
    # pinned-source diagnostic and does not assert a Python parity failure.
    assert first_surge_change > .005
    assert second_surge_change > .03
    assert second_surge_change > 3 * first_surge_change

    assert np.all(work_internal > 0)
    assert np.all(np.abs(np.diff(work_internal)) / work_internal[:-1] < .002)
    # The internal orifice log has strictly increasing solver timestamps;
    # output-grid quadrature remains close enough to audit its work.
    assert np.all(np.abs(work_output - work_internal) / work_internal < .003)


def test_internal_flexible_and_heave_acceleration_integrate_to_speed():
    with np.load(FIXTURES / "owc_orifice_dt_30s/internal_trace.npz",
                 allow_pickle=False) as trace, \
            np.load(FIXTURES / "owc_orifice_dt_30s/source.npz",
                    allow_pickle=False) as source:
        time = trace["time"]
        assert time.ndim == 1 and len(time) > 70_000
        assert time[0] == 0 and time[-1] == 30
        assert np.all(np.diff(time) > 0)
        output_time = source["time"]
        for logged, expected in (
            ("surge_position", source["motion"][0, :, 0]),
            ("heave_velocity", source["velocity"][0, :, 1]),
            ("flex_velocity", source["velocity"][0, :, 3]),
        ):
            np.testing.assert_allclose(
                np.interp(output_time, time, trace[logged]), expected,
                rtol=0, atol=1e-10,
            )

        for speed, acceleration in (
            ("heave_velocity", "heave_acceleration"),
            ("flex_velocity", "flex_hydro_acceleration"),
        ):
            velocity = trace[speed]
            force_acceleration = trace[acceleration]
            integrated = velocity[0] + cumulative_trapezoid(
                force_acceleration, time, initial=0,
            )
            internal_error = float(np.max(np.abs(integrated - velocity)))
            assert internal_error < .02

            # Output-grid quadrature aliases this rapidly varying signal;
            # it cannot be used to infer a missing model force term.
            coarse_velocity = np.interp(output_time, time, velocity)
            coarse_acceleration = np.interp(
                output_time, time, force_acceleration,
            )
            coarse_integrated = (coarse_velocity[0]
                                 + cumulative_trapezoid(
                                     coarse_acceleration, output_time,
                                     initial=0,
                                 ))
            coarse_error = float(np.max(np.abs(coarse_integrated
                                               - coarse_velocity)))
            assert coarse_error > 15 * internal_error


def test_maxstep_only_isolates_most_of_30s_source_surge_shift():
    with np.load(FIXTURES / "owc_orifice_dt_30s/source.npz",
                 allow_pickle=False) as clocks, \
            np.load(FIXTURES / "owc_orifice_dt_30s/maxstep_only.npz",
                    allow_pickle=False) as maxstep, \
            np.load(FIXTURES / "owc_orifice_full_maxstep/source.npz",
                    allow_pickle=False) as full:
        np.testing.assert_array_equal(maxstep["time"], clocks["time"])
        np.testing.assert_array_equal(maxstep["wave"], clocks["wave"][0])
        np.testing.assert_array_equal(maxstep["components"],
                                      clocks["components"][0])
        assert float(maxstep["simu_dt"]) == .005
        assert float(maxstep["max_step"]) == .00125
        assert maxstep["motion"].shape == (6001, 4)
        assert maxstep["velocity"].shape == (6001, 4)

        # This comparison holds the source wave/radiation grid at 5 ms and
        # changes only the initialized ode23t maximum step.
        previous_maxstep = full["half_motion"][:6001, 0]
        surge_maxstep_shift = float(np.max(np.abs(
            maxstep["motion"][:, 0] - previous_maxstep)))
        assert .03 < surge_maxstep_shift < .04

        fine = clocks["motion"][2]
        uniform_half = clocks["motion"][1]
        residual = fine - maxstep["motion"]
        source_gap = fine - uniform_half
        assert np.max(np.abs(residual[:, 0])) < .006
        assert (1 - np.sum(residual[:, 0]**2)
                / np.sum(source_gap[:, 0]**2)) > .98
        # Heave responds differently, so do not generalize surge attribution
        # to all four coordinates.
        assert (1 - np.sum(residual[:, 1]**2)
                / np.sum(source_gap[:, 1]**2)) < .8
        assert abs(float(maxstep["work_internal"])
                   - clocks["work_internal"][2]) \
            / clocks["work_internal"][2] < .001

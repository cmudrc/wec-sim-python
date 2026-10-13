"""Check the opt-in RM3 unilateral PTO stop dynamics."""

import json
from pathlib import Path

import numpy as np
import pytest

from wecsim import LinearHardStops, solve_rm3_regular
from wecsim.dynamics.caseDynamics import run_case


ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / "examples/rm3.json"
HYDRO = ROOT / "examples/data/rm3.h5"
SOURCE_FORCE = (ROOT / "tests/fixtures/matlab"
                / "rm3_end_stops_force_samples.csv")


def test_published_stop_transition_matches_compact_matlab_source():
    samples = np.loadtxt(SOURCE_FORCE, delimiter=",", skiprows=1)
    assert samples.shape == (6, 5)
    np.testing.assert_array_equal(samples[:, 0],
                                  [0, 5166, 8627, 8948, 22725, 28451])
    stroke, speed, source_force = samples[:, 2], samples[:, 3], samples[:, 4]
    stops = LinearHardStops(-0.6, 0.6, 1e8, 1e8)
    predicted = -1_200_000 * speed + stops.force(stroke, speed)
    np.testing.assert_allclose(predicted, source_force, rtol=0, atol=1e-5)

    # Two source samples lie inside the 0.1 mm transition. A sharp spring
    # misses each logged total PTO force by more than 1 kN there.
    transition = (np.abs(stroke) > 0.6) & (np.abs(stroke) < 0.6001)
    assert np.count_nonzero(transition) == 2
    sharp_force = (-1_200_000 * speed
                   + 1e8 * np.maximum(-0.6 - stroke, 0)
                   - 1e8 * np.maximum(stroke - 0.6, 0))
    assert np.min(np.abs(sharp_force[transition]
                         - source_force[transition])) > 1000


def test_contact_damper_only_removes_energy_and_cannot_pull_on_rebound():
    plain = LinearHardStops(-.6, .6, 1000, 1000)
    damped = LinearHardStops(-.6, .6, 1000, 1000,
                              lower_damping=100, upper_damping=100)
    stroke = np.linspace(-.8, .8, 101)[:, None]
    speed = np.linspace(-2, 2, 101)[None, :]
    damping_force = damped.force(stroke, speed) - plain.force(stroke, speed)
    assert np.max(damping_force * speed) <= 1e-12

    # Once the damper exceeds the spring on rebound, the contact opens;
    # neither stop may pull the PTO toward itself.
    assert damped.force(.61, -1) == 0
    assert damped.force(-.61, 1) == 0
    assert damped.force(.61, 1) < 0
    assert damped.force(-.61, -1) > 0


def test_opt_in_stops_are_step_stable_and_report_the_reaction():
    stops = LinearHardStops(-0.6, 0.6, 1e8, 1e8)
    settings = dict(pto_hard_stops=stops, end_time=10, ramp_time=0)
    coarse = solve_rm3_regular(HYDRO, dt=0.1, **settings)
    fine = solve_rm3_regular(HYDRO, dt=0.05, **settings)
    assert np.min(coarse.pto_stroke) < -0.6
    assert np.max(coarse.pto_stroke) > 0.6
    assert np.max(np.abs(coarse.pto_stroke - fine.pto_stroke[::2])) < 5e-5
    np.testing.assert_allclose(
        coarse.pto_stop_force,
        stops.force(coarse.pto_stroke, coarse.pto_velocity),
        rtol=0, atol=1e-9,
    )
    np.testing.assert_allclose(
        coarse.pto_force,
        -1_200_000 * coarse.pto_velocity + coarse.pto_stop_force,
        rtol=0, atol=1e-9,
    )


def test_case_runner_accepts_explicit_stop_settings():
    case = json.loads(EXAMPLE.read_text(encoding="utf-8"))
    case["simulation"].update(end_time=10, ramp_time=0)
    case["pto"]["hard_stops"] = {
        "lower_bound": -0.6, "upper_bound": 0.6,
        "lower_stiffness": 1e8, "upper_stiffness": 1e8,
    }
    response = run_case(case, base_dir=EXAMPLE.parent)
    extras = dict(response.extra_outputs)
    assert response.pto_force.shape == (101,)
    assert extras["pto_stop_force"].shape == (101,)
    assert np.isfinite(response.body_position).all()
    assert np.max(np.abs(extras["pto_stop_force"])) > 1e6


def test_stops_reject_the_source_delay_scheme():
    stops = LinearHardStops(-0.6, 0.6, 1e8, 1e8)
    with pytest.raises(ValueError, match="implicit added mass"):
        solve_rm3_regular(
            HYDRO, pto_hard_stops=stops,
            radiation_memory=2, added_mass_scheme="simulink_delay",
            end_time=0,
        )

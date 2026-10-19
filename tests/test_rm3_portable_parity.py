"""Keep a published MATLAB RM3 trajectory paired in ordinary Python CI."""

from pathlib import Path

import numpy as np
import pytest

from wecsim import RegularWave, WEC
from wecsim.dynamics.linearHeave import solve_two_body_regular_heave


ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT / "tests/fixtures/matlab/rm3_regular_wave.npz"
HYDRO = ROOT / "examples/data/rm3.h5"
WAVE_MARKERS_HYDRO = ROOT / "tests/fixtures/matlab/rm3_mooring_matrix/rm3.h5"


def _max_error(actual, expected, limit, label):
    assert actual.shape == expected.shape, f"{label}: mismatched shape"
    assert np.isfinite(actual).all() and np.isfinite(expected).all(), label
    error = float(np.max(np.abs(actual - expected)))
    assert error < limit, f"{label}: {error:.6g} exceeds {limit:.6g}"


@pytest.mark.parametrize(
    "hydro", [HYDRO, WAVE_MARKERS_HYDRO],
    ids=["bundled-core", "published-wave-markers"],
)
def test_published_rm3_regular_wave_body_and_pto_motion(hydro):
    """Pair both the bundled and published-input hydrodynamics."""
    wec = WEC("RM3 published regular wave")
    float_body = wec.body("float", hydro, inertia=(0, 21_306_090.66, 0))
    spar = wec.body("spar", hydro, inertia=(0, 94_407_091.24, 0))
    wec.floating_joint(float_body, spar, damping=1_200_000,
                       radiation_method="constant")
    result = wec.run(RegularWave(2.5, 8), dt=0.1, end_time=400,
                     ramp_time=100)

    with np.load(REFERENCE, allow_pickle=False) as source:
        assert source["time"].shape == (4001,)
        assert source["body_position"].shape == (2, 4001, 6)
        assert source["body_velocity"].shape == (2, 4001, 6)
        np.testing.assert_allclose(result.time, source["time"], rtol=0,
                                   atol=1e-10)
        body_position = source["body_position"]
        body_velocity = source["body_velocity"]

        # Source columns follow [surge, sway, heave, roll, pitch, yaw].
        # The published joint has no sway, roll, or yaw motion.
        position_limits = (0.04, 1e-10, 0.006, 1e-10, 0.001, 1e-10)
        velocity_limits = (0.03, 1e-10, 0.004, 1e-10, 0.001, 1e-10)
        for body_index, name in enumerate(("float", "spar")):
            for dof, (position_limit, velocity_limit) in enumerate(
                zip(position_limits, velocity_limits)
            ):
                motion = result.bodies[name]
                _max_error(motion.position[:, dof],
                           body_position[body_index, :, dof],
                           position_limit, f"{name} position DOF {dof}")
                _max_error(motion.velocity[:, dof],
                           body_velocity[body_index, :, dof],
                           velocity_limit, f"{name} velocity DOF {dof}")

        pto = result.ptos["relative_heave"]
        for label, actual, expected, limit in (
            ("stroke", pto.stroke, source["pto_stroke"], 0.006),
            ("speed", pto.velocity, source["pto_speed"], 0.005),
            ("force", pto.force, source["pto_force"], 6_000),
            # MATLAB logs mechanical power with the opposite sign from
            # Python's absorbed-power convention for this linear damper.
            ("absorbed power", pto.absorbed_power, -source["pto_power"], 7_000),
        ):
            _max_error(actual, expected, limit, f"RM3 PTO {label}")

        # Pointwise power gates do not bound an accumulated output over 400 s.
        # Check the complete run and separate ramp/steady intervals so errors
        # in one part cannot cancel those in another.
        source_absorbed_power = -source["pto_power"]
        for start, stop in ((0, 400), (0, 100), (100, 250), (250, 400)):
            window = (result.time >= start) & (result.time <= stop)
            source_work = np.trapezoid(
                source_absorbed_power[window], result.time[window],
            )
            python_work = np.trapezoid(
                pto.absorbed_power[window], result.time[window],
            )
            assert source_work > 1e6
            relative_gap = abs(python_work - source_work) / source_work
            assert relative_gap < .01, (
                f"RM3 PTO absorbed work {start}–{stop} s: "
                f"{relative_gap:.4%} exceeds 1%"
            )


def test_published_rm3_heave_excitation_and_relative_pto():
    """Pair the reduced heave solver's independent wave and PTO force paths."""
    response = solve_two_body_regular_heave(
        HYDRO, wave_height=2.5, wave_period=8,
        pto_damping=1_200_000,
    )
    with np.load(REFERENCE, allow_pickle=False) as source:
        np.testing.assert_allclose(response.time, source["time"], rtol=0,
                                   atol=1e-10)
        expected_excitation = source["body_heave_excitation_force"]
        assert expected_excitation.shape == (2, 4001)
        for body in range(2):
            _max_error(response.position[:, body],
                       source["body_position"][body, :, 2], 0.0065,
                       f"heave body {body + 1} position")
            _max_error(response.velocity[:, body],
                       source["body_velocity"][body, :, 2], 0.0065,
                       f"heave body {body + 1} speed")
            np.testing.assert_allclose(
                response.excitation_force[:, body], expected_excitation[body],
                rtol=0, atol=1e-5,
            )

        stroke = response.displacement[:, 0] - response.displacement[:, 1]
        speed = response.velocity[:, 0] - response.velocity[:, 1]
        _max_error(stroke, source["pto_stroke"], 0.009, "heave PTO stroke")
        _max_error(speed, source["pto_speed"], 0.0075, "heave PTO speed")
        _max_error(response.pto_force, source["pto_force"], 9_000,
                   "heave PTO force")
        _max_error(response.pto_force * speed, source["pto_power"], 12_000,
                   "heave PTO mechanical power")

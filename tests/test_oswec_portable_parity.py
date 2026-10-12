"""Keep the published directional OSWEC run paired in ordinary Python CI."""

from pathlib import Path

import numpy as np
import pytest

from wecsim import PMWave, WEC, WorldPoint
from wecsim.environment.irregularWave import (
    pm_equal_energy_components,
    synthesize_irregular_response,
)


ROOT = Path(__file__).resolve().parents[1]
HYDRO = ROOT / "examples/data/oswec.h5"
REFERENCE = ROOT / "tests/fixtures/matlab/oswec_irregular_wave.npz"


@pytest.fixture(scope="module")
def source():
    with np.load(REFERENCE, allow_pickle=False) as saved:
        arrays = {name: saved[name] for name in saved.files}
    assert arrays["time"].shape == (4001,)
    assert arrays["wave_components"].shape == (500, 6)
    assert arrays["wave_directions"].shape == (3, 2)
    assert arrays["body_position"].shape == (2, 4001, 6)
    assert arrays["body_velocity"].shape == (2, 4001, 6)
    assert arrays["flap_excitation_force"].shape == (4001, 6)
    return arrays


def _max_error(actual, expected, limit, label):
    assert actual.shape == expected.shape, f"{label}: mismatched shape"
    assert np.isfinite(actual).all() and np.isfinite(expected).all(), label
    error = float(np.max(np.abs(actual - expected)))
    assert error < limit, f"{label}: {error:.6g} exceeds {limit:.6g}"


def test_published_oswec_irregular_wave_and_excitation(source):
    components = pm_equal_energy_components(
        HYDRO, significant_height=2.5, peak_period=8,
        directions=source["wave_directions"][:, 0],
        spreading=source["wave_directions"][:, 1],
        phase=source["wave_components"][:, 3:],
    )
    for actual, expected in (
        (components.omega, source["wave_components"][:, 0]),
        (components.spectral_amplitude, source["wave_components"][:, 1]),
        (components.d_omega, source["wave_components"][:, 2]),
    ):
        np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-12)
    response = synthesize_irregular_response(
        HYDRO, components, dt=.1, end_time=400, ramp_time=100,
    )
    np.testing.assert_allclose(response.time, source["time"], rtol=0,
                               atol=1e-10)
    _max_error(response.elevation, source["wave_elevation"], 1e-11,
               "OSWEC wave elevation")
    _max_error(response.excitation_force, source["flap_excitation_force"],
               1e-6, "OSWEC flap excitation")


def test_public_oswec_irregular_motion_and_pto(source, tmp_path):
    phase_file = tmp_path / "phase.csv"
    # The published input has no phase seed. Replay its recorded realization.
    np.savetxt(phase_file, source["wave_components"][:, 3:], delimiter=",")
    wec = WEC("Published directional OSWEC")
    flap = wec.body("flap", HYDRO, mass=127_000, inertia=(1.85e6,) * 3)
    base = wec.body("base", HYDRO, mass=999, inertia=(999,) * 3)
    wec.fixed_hinge(flap, base, location=WorldPoint(0, 0, -10),
                    pto_location=WorldPoint(0, 0, -8.9), damping=12_000)
    result = wec.run(
        PMWave(2.5, 8, directions=tuple(source["wave_directions"][:, 0]),
               spreading=tuple(source["wave_directions"][:, 1]),
               phase_file=phase_file),
        dt=.1, end_time=400, ramp_time=100, radiation_memory=30,
    )
    np.testing.assert_allclose(result.time, source["time"], rtol=0,
                               atol=1e-10)
    _max_error(result.wave_elevation, source["wave_elevation"], 1e-11,
               "public OSWEC wave elevation")
    position_limits = (.02, 1e-10, .005, 1e-10, .004, 1e-10)
    velocity_limits = (.03, 1e-10, .005, 1e-10, .005, 1e-10)
    for body_index, name in enumerate(("flap", "base")):
        motion = result.bodies[name]
        for dof in range(6):
            position_limit = position_limits[dof] if body_index == 0 else 1e-10
            velocity_limit = velocity_limits[dof] if body_index == 0 else 1e-10
            _max_error(motion.position[:, dof],
                       source["body_position"][body_index, :, dof],
                       position_limit, f"OSWEC {name} position DOF {dof}")
            _max_error(motion.velocity[:, dof],
                       source["body_velocity"][body_index, :, dof],
                       velocity_limit, f"OSWEC {name} speed DOF {dof}")
    pto = result.ptos["hinge"]
    for label, actual, expected, limit in (
        ("stroke", pto.stroke, source["pto_stroke"], .004),
        ("speed", pto.velocity, source["pto_speed"], .005),
        ("torque", pto.force, source["pto_force"], 60),
        # MATLAB reports mechanical power with the opposite sign.
        ("power", -pto.absorbed_power, source["pto_power"],
         max(25, .06 * np.max(np.abs(source["pto_power"])))),
    ):
        _max_error(actual, expected, limit, f"OSWEC PTO {label}")

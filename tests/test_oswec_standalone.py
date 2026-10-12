"""Exercise the published OSWEC settings without any MATLAB-generated output."""

from pathlib import Path

import numpy as np
import pytest

from wecsim.dynamics.hingePitch import solve_hinged_pitch_from_excitation
from wecsim.dynamics.caseDynamics import run_case
from wecsim.environment.irregularWave import (
    pm_equal_energy_components, synthesize_irregular_response,
)


def test_seeded_oswec_wave_and_pitch_pipeline():
    h5_file = Path(__file__).parent / "test_objects" / "test_bodyclass" / "testData" / "hydroData" / "oswec.h5"
    settings = dict(
        significant_height=2.5, peak_period=8,
        directions=[0, 30, 90], spreading=[0.1, 0.2, 0.7], seed=7,
    )
    components = pm_equal_energy_components(h5_file, **settings)
    replay = pm_equal_energy_components(h5_file, **settings)
    np.testing.assert_array_equal(components.phase, replay.phase)
    wave = synthesize_irregular_response(
        h5_file, components, dt=0.1, end_time=40, ramp_time=10,
    )
    motion = solve_hinged_pitch_from_excitation(
        h5_file, wave.excitation_force, hinge_z=-8.9,
        body_mass=127_000, pitch_inertia=1.85e6, pto_damping=12_000,
    )
    assert len(motion.time) == 401
    assert np.isfinite(wave.elevation).all()
    assert np.isfinite(wave.excitation_force).all()
    assert np.isfinite(motion.angle).all()
    assert 0 < np.max(np.abs(motion.angle)) < 1
    np.testing.assert_allclose(motion.time, wave.time, rtol=0, atol=1e-12)


def test_regular_hinged_flap_reports_fixed_nonhydrodynamic_base():
    h5_file = (Path(__file__).parent / "test_objects" / "test_bodyclass"
               / "testData" / "hydroData" / "oswec.h5")
    case = {
        "simulation": {"dt": 0.1, "end_time": 8, "ramp_time": 0},
        "wave": {"type": "regular", "height": 2.5, "period": 8},
        "bodies": [
            {"hydro_file": str(h5_file.resolve()), "mass": 127_000,
             "pitch_inertia": 1.85e6},
            {"name": "base", "nonhydro": True, "fixed": True,
             "center_gravity": [0, 0, -10.9], "mass": 999,
             "inertia": [1, 1, 1]},
        ],
        "constraint": {"kind": "fixed_hinge", "location": [0, 0, -10]},
        "pto": {"kind": "pitch", "damping": 0,
                "location": [0, 0, -8.9]},
    }
    response = run_case(case)
    assert response.body_position.shape == (81, 2, 6)
    assert len(response.hydro_files) == 1
    np.testing.assert_array_equal(
        response.body_position[:, 1, :],
        np.tile([0, 0, -10.9, 0, 0, 0], (81, 1)),
    )
    assert np.max(np.abs(response.body_position[:, 0, 4])) > 0.01
    np.testing.assert_array_equal(response.pto_force, np.zeros(81))
    flap_only = run_case({**case, "bodies": case["bodies"][:1]})
    np.testing.assert_array_equal(
        response.body_position[:, 0, :], flap_only.body_position[:, 0, :],
    )
    with pytest.raises(ValueError, match="needs pto.location"):
        run_case({**case, "pto": {"kind": "pitch", "damping": 0}})

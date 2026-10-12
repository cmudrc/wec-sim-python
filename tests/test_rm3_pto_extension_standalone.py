"""Exercise no-wave RM3 floating-joint initial states without MATLAB."""

from pathlib import Path

import numpy as np
import pytest

from wecsim.dynamics.caseDynamics import run_case


HYDRO = (Path(__file__).resolve().parents[1] / "examples/data/rm3.h5")


def _case():
    return {
        "simulation": {"dt": 0.1, "end_time": 3,
                       "radiation_memory": 15},
        "wave": {"type": "none"},
        "bodies": [
            {"hydro_file": str(HYDRO), "hydro_body": 1,
             "mass": "equilibrium", "pitch_inertia": 21_306_090.66},
            {"hydro_file": str(HYDRO), "hydro_body": 2,
             "mass": "equilibrium", "pitch_inertia": 94_407_091.24},
        ],
        "constraint": {
            "kind": "floating_joint", "location": [0, 0, 0],
            "initial_coordinate": {"float_heave": 5},
        },
        "pto": {"kind": "relative_heave", "damping": 0},
    }


def test_no_wave_initial_float_displacement_drives_free_decay():
    case = _case()
    result = run_case(case)
    equilibrium = run_case({
        **case,
        "simulation": {**case["simulation"], "end_time": 0},
        "constraint": {**case["constraint"], "initial_coordinate": {}},
    })
    assert result.body_position.shape == (31, 2, 6)
    assert result.wave_elevation is None
    np.testing.assert_allclose(
        result.body_position[0, 0, 2] - equilibrium.body_position[0, 0, 2],
        5, rtol=0, atol=1e-10,
    )
    np.testing.assert_allclose(
        result.body_position[0, 1, 2], equilibrium.body_position[0, 1, 2],
        rtol=0, atol=1e-10,
    )
    assert np.max(np.abs(result.body_velocity[1:, :, 2])) > 0.01
    assert np.isfinite(result.body_position).all()
    assert np.isfinite(result.body_velocity).all()


def test_initial_coordinate_validation():
    case = _case()
    case["constraint"]["initial_coordinate"] = {"unsupported": 5}
    with pytest.raises(ValueError, match="unknown coordinates"):
        run_case(case)
    case["constraint"]["initial_coordinate"] = [0, 5]
    with pytest.raises(ValueError, match="one finite value per coordinate"):
        run_case(case)


def test_named_initial_speed_reaches_both_body_velocities():
    case = _case()
    case["simulation"]["end_time"] = 0
    case["constraint"]["initial_coordinate"] = {}
    case["constraint"]["initial_speed"] = {
        "surge": 0.2, "float_heave": -0.3, "spar_heave": 0.4,
    }
    result = run_case(case)
    np.testing.assert_allclose(result.body_velocity[0, 0, [0, 2]],
                               [0.2, -0.3], atol=1e-12)
    np.testing.assert_allclose(result.body_velocity[0, 1, [0, 2]],
                               [0.2, 0.4], atol=1e-12)

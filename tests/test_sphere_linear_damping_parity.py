"""Pair a body damping setting with a locally generated WEC-Sim Sphere run."""

from pathlib import Path

import numpy as np
import pytest

from wecsim import NoWave, WEC


ROOT = Path(__file__).resolve().parents[1]
HYDRO = ROOT / "examples/data/sphere.h5"
SOURCE = ROOT / "tests/fixtures/matlab/sphere_linear_damping/source.npz"
PUBLISHED = ROOT / "tests/fixtures/matlab/sphere_free_decay.npz"


def test_body_linear_damping_matches_derived_sphere_free_decay():
    with np.load(SOURCE, allow_pickle=False) as saved:
        time = saved["time"].copy()
        position = saved["position"].copy()
        velocity = saved["velocity"].copy()
        damping = saved["linear_damping"].copy()
    with np.load(PUBLISHED, allow_pickle=False) as published:
        np.testing.assert_allclose(time, published["time"], rtol=0, atol=1e-12)
        # The added body damper changes the source motion substantially.
        assert np.max(np.abs(position[:, 2] - published["body_position"][1, :, 2])) > .3

    wec = WEC("Sphere with body heave damping")
    sphere = wec.body("sphere", HYDRO, linear_damping=damping)
    wec.coordinate("heave", sphere.move("heave"))
    result = wec.run(NoWave(), dt=.01, end_time=40,
                     radiation_memory=15, initial_coordinate={"heave": 1})
    np.testing.assert_allclose(result.time, time, rtol=0, atol=1e-12)
    assert result.case["bodies"][0]["linear_damping"] == damping.tolist()
    actual = result.bodies["sphere"]
    for dof in range(6):
        position_limit = 1.2e-4 if dof == 2 else 1e-10
        speed_limit = 1.6e-4 if dof == 2 else 1e-10
        assert np.max(np.abs(actual.position[:, dof] - position[:, dof])) < position_limit
        assert np.max(np.abs(actual.velocity[:, dof] - velocity[:, dof])) < speed_limit


@pytest.mark.parametrize("matrix", [np.ones(6), np.full((6, 6), np.nan)])
def test_body_rejects_invalid_damping_without_reserving_name(matrix):
    wec = WEC("body damping configuration")
    with pytest.raises(ValueError, match="linear_damping needs a finite 6x6 matrix"):
        wec.body("sphere", HYDRO, linear_damping=matrix)
    assert wec.bodies == []
    wec.body("sphere", HYDRO, linear_damping=np.zeros((6, 6)))
    assert len(wec.bodies) == 1


def test_specialized_joint_rejects_unsupported_body_damping():
    wec = WEC("RM3 with custom damping")
    hydro = ROOT / "examples/data/rm3.h5"
    float_body = wec.body("float", hydro, inertia=(0, 21_306_090.66, 0),
                          linear_damping=np.eye(6))
    spar = wec.body("spar", hydro, inertia=(0, 94_407_091.24, 0))
    wec.floating_joint(float_body, spar)
    with pytest.raises(ValueError, match="without extra force models"):
        wec.to_case(NoWave(), dt=.1, end_time=.2)

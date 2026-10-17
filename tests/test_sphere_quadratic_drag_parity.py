"""Pair public body quadratic drag with a local WEC-Sim Sphere source run."""

from pathlib import Path

import numpy as np
import pytest

from wecsim import NoWave, WEC, WorldPoint


ROOT = Path(__file__).resolve().parents[1]
HYDRO = ROOT / "examples/data/sphere.h5"
SOURCE = ROOT / "tests/fixtures/matlab/sphere_quadratic_drag/source.npz"
SOURCE_CROSS = ROOT / "tests/fixtures/matlab/sphere_quadratic_drag/source_cross.npz"
PUBLISHED = ROOT / "tests/fixtures/matlab/sphere_free_decay.npz"


def test_body_quadratic_drag_matches_derived_sphere_free_decay():
    with np.load(SOURCE, allow_pickle=False) as saved:
        time = saved["time"].copy()
        position = saved["position"].copy()
        velocity = saved["velocity"].copy()
        force = saved["viscous_force"].copy()
        drag = saved["quadratic_drag"].copy()
    with np.load(PUBLISHED, allow_pickle=False) as published:
        np.testing.assert_allclose(time, published["time"], rtol=0, atol=1e-12)
        # This is a distinct damped trajectory, not an all-zero setting.
        assert np.max(np.abs(position[:, 2] - published["body_position"][1, :, 2])) > .4

    # WEC-Sim logs the positive damping term; its body block subtracts it
    # from applied force. This also checks the matrix-times-v|v| convention.
    source_drag = (drag @ (velocity * np.abs(velocity)).T).T
    assert np.max(np.abs(source_drag - force)) < 1e-5

    wec = WEC("Sphere with quadratic body heave drag")
    sphere = wec.body("sphere", HYDRO, quadratic_drag=drag)
    wec.coordinate("heave", sphere.move("heave"))
    result = wec.run(NoWave(), dt=.01, end_time=40,
                     radiation_memory=15, initial_coordinate={"heave": 1})
    np.testing.assert_allclose(result.time, time, rtol=0, atol=1e-12)
    assert result.case["bodies"][0]["quadratic_drag"] == drag.tolist()
    actual = result.bodies["sphere"]
    for dof in range(6):
        position_limit = 1e-4 if dof == 2 else 1e-10
        speed_limit = 1e-4 if dof == 2 else 1e-10
        assert np.max(np.abs(actual.position[:, dof] - position[:, dof])) < position_limit
        assert np.max(np.abs(actual.velocity[:, dof] - velocity[:, dof])) < speed_limit


def test_off_diagonal_quadratic_drag_matches_coupled_sphere_motion():
    with np.load(SOURCE_CROSS, allow_pickle=False) as saved:
        time = saved["time"].copy()
        position = saved["position"].copy()
        velocity = saved["velocity"].copy()
        force = saved["viscous_force"].copy()
        drag = saved["quadratic_drag"].copy()
    with np.load(SOURCE, allow_pickle=False) as diagonal:
        # A heave-speed term in the surge-force row exercises the cross entry.
        assert np.max(np.abs(position[:, 0] - diagonal["position"][:, 0])) > 4
    source_drag = (drag @ (velocity * np.abs(velocity)).T).T
    assert np.max(np.abs(source_drag - force)) < 1e-5

    wec = WEC("Sphere with cross-coordinate body drag")
    sphere = wec.body(
        "sphere", HYDRO, inertia=(20_907_301, 21_306_090.66, 37_085_481.11),
        quadratic_drag=drag,
    )
    wec.coordinate("surge", sphere.move("surge"))
    wec.coordinate("heave", sphere.move("heave"))
    wec.coordinate("pitch", sphere.move("pitch", pivot=WorldPoint(0, 0, 0)))
    result = wec.run(NoWave(), dt=.01, end_time=40,
                     radiation_memory=15, initial_coordinate={"heave": 1})
    np.testing.assert_allclose(result.time, time, rtol=0, atol=1e-12)
    assert result.case["bodies"][0]["quadratic_drag"] == drag.tolist()
    actual = result.bodies["sphere"]
    for dof, position_limit, speed_limit in (
        (0, 2e-4, 2e-5), (2, 1e-4, 1e-4), (4, 3e-7, 3e-7),
        (1, 1e-10, 1e-10), (3, 1e-10, 1e-10), (5, 1e-10, 1e-10),
    ):
        assert np.max(np.abs(actual.position[:, dof] - position[:, dof])) < position_limit
        assert np.max(np.abs(actual.velocity[:, dof] - velocity[:, dof])) < speed_limit


@pytest.mark.parametrize("matrix", [np.ones(6), np.full((6, 6), np.nan)])
def test_body_rejects_invalid_quadratic_drag_without_reserving_name(matrix):
    wec = WEC("body quadratic drag configuration")
    with pytest.raises(ValueError, match="quadratic_drag needs a finite 6x6 matrix"):
        wec.body("sphere", HYDRO, quadratic_drag=matrix)
    assert wec.bodies == []
    wec.body("sphere", HYDRO, quadratic_drag=np.zeros((6, 6)))
    assert len(wec.bodies) == 1


def test_specialized_joint_rejects_unsupported_body_quadratic_drag():
    wec = WEC("RM3 with custom quadratic drag")
    hydro = ROOT / "examples/data/rm3.h5"
    float_body = wec.body("float", hydro, inertia=(0, 21_306_090.66, 0),
                          quadratic_drag=np.eye(6))
    spar = wec.body("spar", hydro, inertia=(0, 94_407_091.24, 0))
    wec.floating_joint(float_body, spar)
    with pytest.raises(ValueError, match="without extra force models"):
        wec.to_case(NoWave(), dt=.1, end_time=.2)

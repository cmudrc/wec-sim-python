"""Pair a public hydrostatic override with an independent MATLAB Sphere run."""

from pathlib import Path

import numpy as np
import pytest

from wecsim import NoWave, WEC


ROOT = Path(__file__).resolve().parents[1]
HYDRO = ROOT / "examples/data/sphere.h5"
SOURCE = ROOT / "tests/fixtures/matlab/sphere_hydrostatic_stiffness/source.npz"
HALF_STEP = ROOT / "tests/fixtures/matlab/sphere_hydrostatic_stiffness/half_step.npz"
PUBLISHED = ROOT / "tests/fixtures/matlab/sphere_free_decay.npz"


def test_body_hydrostatic_stiffness_matches_derived_sphere_free_decay():
    with np.load(SOURCE, allow_pickle=False) as saved:
        time = saved["time"].copy()
        position = saved["position"].copy()
        velocity = saved["velocity"].copy()
        stiffness = saved["hydrostatic_stiffness"].copy()
    with np.load(PUBLISHED, allow_pickle=False) as published:
        np.testing.assert_allclose(time, published["time"], rtol=0, atol=1e-12)
        assert np.max(np.abs(position[:, 2] - published["body_position"][1, :, 2])) > .1

    wec = WEC("Sphere with heave restoring override")
    sphere = wec.body("sphere", HYDRO, hydrostatic_stiffness=stiffness)
    wec.coordinate("heave", sphere.move("heave"))
    result = wec.run(NoWave(), dt=.01, end_time=40,
                     radiation_memory=15, initial_coordinate={"heave": 1},
                     solver="wecsim_ode4")
    np.testing.assert_allclose(result.time, time, rtol=0, atol=1e-12)
    assert result.case["bodies"][0]["hydrostatic_stiffness"] == stiffness.tolist()
    actual = result.bodies["sphere"]
    for dof in range(6):
        position_limit = 1e-8 if dof == 2 else 1e-10
        speed_limit = 1e-8 if dof == 2 else 1e-10
        assert np.max(np.abs(actual.position[:, dof] - position[:, dof])) < position_limit
        assert np.max(np.abs(actual.velocity[:, dof] - velocity[:, dof])) < speed_limit


def test_physical_solver_pairs_with_matched_half_step_sphere_source():
    """Resolve the physical solver's published-step miss without changing its force law."""
    with np.load(SOURCE, allow_pickle=False) as coarse:
        coarse_time = coarse["time"].copy()
        coarse_position = coarse["position"].copy()
        coarse_velocity = coarse["velocity"].copy()
        stiffness = coarse["hydrostatic_stiffness"].copy()
    with np.load(HALF_STEP, allow_pickle=False) as refined:
        time = refined["time"].copy()
        position = refined["position"].copy()
        velocity = refined["velocity"].copy()
    assert time.shape == (8001,) and position.shape == velocity.shape == (8001, 6)
    np.testing.assert_allclose(time, np.arange(8001) * .005, rtol=0, atol=1e-10)
    np.testing.assert_allclose(time[::2], coarse_time, rtol=0, atol=1e-10)
    assert np.isfinite(position).all() and np.isfinite(velocity).all()

    # The MATLAB source itself changes less than the original 0.12 mm and
    # 0.16 mm/s motion limits when only its 0.01 s step is halved.
    assert np.max(np.abs(position[::2, 2] - coarse_position[:, 2])) < 0.00012
    assert np.max(np.abs(velocity[::2, 2] - coarse_velocity[:, 2])) < 0.00016

    wec = WEC("Sphere matched-step restoring audit")
    sphere = wec.body("sphere", HYDRO, hydrostatic_stiffness=stiffness)
    wec.coordinate("heave", sphere.move("heave"))
    result = wec.run(NoWave(), dt=.005, end_time=40,
                     radiation_memory=15, initial_coordinate={"heave": 1})
    np.testing.assert_allclose(result.time, time, rtol=0, atol=1e-10)
    actual = result.bodies["sphere"]
    for dof in range(6):
        position_limit = 0.00012 if dof == 2 else 1e-10
        speed_limit = 0.00016 if dof == 2 else 1e-10
        assert np.max(np.abs(actual.position[:, dof] - position[:, dof])) < position_limit
        assert np.max(np.abs(actual.velocity[:, dof] - velocity[:, dof])) < speed_limit


@pytest.mark.parametrize("matrix", [np.ones(6), np.full((6, 6), np.nan)])
def test_body_rejects_invalid_hydrostatic_stiffness(matrix):
    wec = WEC("hydrostatic override validation")
    with pytest.raises(ValueError, match="hydrostatic_stiffness needs a finite 6x6 matrix"):
        wec.body("sphere", HYDRO, hydrostatic_stiffness=matrix)
    assert wec.bodies == []


def test_zero_hydrostatic_matrix_uses_hdf5_restoring():
    results = []
    for stiffness in (None, np.zeros((6, 6))):
        wec = WEC("HDF5 restoring fallback")
        sphere = wec.body("sphere", HYDRO, hydrostatic_stiffness=stiffness)
        wec.coordinate("heave", sphere.move("heave"))
        results.append(wec.run(
            NoWave(), dt=.01, end_time=.2, radiation_memory=15,
            initial_coordinate={"heave": 1},
        ))
    np.testing.assert_array_equal(
        results[0].bodies["sphere"].position,
        results[1].bodies["sphere"].position,
    )
    np.testing.assert_array_equal(
        results[0].bodies["sphere"].velocity,
        results[1].bodies["sphere"].velocity,
    )


def test_specialized_joint_rejects_unsupported_hydrostatic_stiffness():
    wec = WEC("RM3 with custom restoring")
    hydro = ROOT / "examples/data/rm3.h5"
    float_body = wec.body("float", hydro, inertia=(0, 21_306_090.66, 0),
                          hydrostatic_stiffness=np.eye(6))
    spar = wec.body("spar", hydro, inertia=(0, 94_407_091.24, 0))
    wec.floating_joint(float_body, spar)
    with pytest.raises(ValueError, match="without extra force models"):
        wec.to_case(NoWave(), dt=.1, end_time=.2)

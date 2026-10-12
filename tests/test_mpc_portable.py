"""Short Sphere MPC checks using the bundled hydro file and a generated fit."""

from pathlib import Path

import numpy as np
import pytest
from scipy.io import savemat

from wecsim.control.mpcControl import (
    integrate_sphere_mpc_force, predict_sphere_excitation,
    simulate_sphere_mpc, solve_sphere_mpc_qp,
)
from wecsim.control.mpcPlant import (
    build_sphere_mpc_matrices, replay_sphere_mpc_plant,
)
from wecsim.control.sphereMpc import run_sphere_mpc


SPHERE = Path(__file__).resolve().parents[1] / "examples/data/sphere.h5"


def _matrices(tmp_path):
    coefficients = tmp_path / "radiation.mat"
    # A stable fourth-order radiation fit gives a small, independent MPC
    # plant without relying on the published MATLAB coefficient artifact.
    savemat(coefficients, {"coeff": {
        "KradNum": [0, 0, 0, 0], "KradDen": [1, 4, 6, 4, 1],
    }})
    return build_sphere_mpc_matrices(
        SPHERE, coefficients, prediction_step=.1, prediction_horizon=.3,
    )


def test_mpc_forecast_follows_a_geometric_force_history():
    history = 100 * .9 ** np.arange(9)
    predicted = predict_sphere_excitation(
        history, history_steps=6, order=2, horizon_steps=3,
    )
    np.testing.assert_allclose(predicted, history[-1] * .9 ** np.arange(4),
                               rtol=0, atol=1e-10)
    with pytest.raises(ValueError, match="finite ordered"):
        predict_sphere_excitation([1, 2, np.nan])


def test_mpc_prediction_matrices_and_plant_replay_obey_force_integrator(tmp_path):
    matrices = _matrices(tmp_path)
    assert matrices.Sx.shape == (9, 9)
    assert matrices.Su.shape == matrices.Sv.shape == (9, 4)
    assert np.linalg.eigvalsh(matrices.H).min() > 0
    np.testing.assert_allclose(matrices.C @ matrices.Bu, [0, 0, 1])
    force = np.array([0., 1., 2., 3.])
    replay = replay_sphere_mpc_plant(
        matrices, force, np.zeros_like(force), dt=.1,
    )
    np.testing.assert_allclose(replay[:, 2], force, rtol=0, atol=1e-12)
    assert np.isfinite(replay).all()
    with pytest.raises(ValueError, match="equal finite force"):
        replay_sphere_mpc_plant(matrices, force, np.zeros(3), dt=.1)


def test_mpc_rate_plan_obeys_limits_and_delayed_force_application(tmp_path):
    matrices = _matrices(tmp_path)
    plan, feasible = solve_sphere_mpc_qp(
        matrices, np.zeros(9), np.full(4, 1e4),
        max_force_rate=5e4,
    )
    assert feasible and np.isfinite(plan).all()
    assert np.max(np.abs(plan)) <= 5e4 + 1e-6
    assert np.linalg.norm(plan) > 0
    with pytest.raises(ValueError, match="positive limits"):
        solve_sphere_mpc_qp(matrices, np.zeros(9), np.zeros(4),
                            max_force_rate=0)

    history = integrate_sphere_mpc_force(
        np.array([10., 20., 30., 40.]), dt=.1, end_time=.5,
        control_step=.2, command_delay=.1,
    )
    np.testing.assert_allclose(history, [0, 0, 1, 2, 4, 6])
    with pytest.raises(ValueError, match="align"):
        integrate_sphere_mpc_force(np.ones(3), dt=.1, end_time=.5,
                                   control_step=.25)


def test_mpc_zero_excitation_remains_at_rest_through_active_controller(tmp_path):
    matrices = _matrices(tmp_path)
    response = simulate_sphere_mpc(
        matrices, np.zeros(11), dt=.1, control_step=.1,
        command_delay=.1, start_time=.5, history_steps=3, order=1,
    )
    assert response.feasible.all()
    np.testing.assert_array_equal(response.internal_state, 0)
    np.testing.assert_array_equal(response.pto_force, 0)
    np.testing.assert_array_equal(response.command_rate, 0)


def test_short_public_sphere_mpc_reports_signed_pto_power(tmp_path):
    coefficients = tmp_path / "radiation.mat"
    savemat(coefficients, {"coeff": {
        "KradNum": [0, 0, 0, 0], "KradDen": [1, 4, 6, 4, 1],
    }})
    result = run_sphere_mpc(
        SPHERE, coefficients, significant_height=.5,
        frequency_count=8, seed=4, dt=.01, end_time=.1, ramp_time=0,
        radiation_memory=.05, control_step=.01, prediction_horizon=.03,
        command_delay=.01, start_time=.05, history_steps=3, order=1,
    )
    assert len(result.time) == 11
    assert np.isfinite(result.position).all()
    assert result.feasible.all()
    np.testing.assert_allclose(result.absorbed_power,
                               -result.pto_force * result.velocity)
    np.testing.assert_array_equal(result.pto_force[:6], 0)

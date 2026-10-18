"""Compare the published Sphere MPC optimizer on Python-generated inputs."""

from pathlib import Path

import numpy as np

from wecsim.control.mpcControl import solve_sphere_mpc_qp
from wecsim.control.mpcPlant import build_sphere_mpc_matrices


ROOT = Path(__file__).parent / "fixtures/matlab"
MODEL = ROOT / "sphere_mpc"
FIXTURE = ROOT / "sphere_mpc_independent_qp/source.npz"


def test_matlab_quadprog_on_independent_python_states_and_forecasts():
    """Gate all active commands, including the force-rate saturation."""
    with np.load(FIXTURE, allow_pickle=False) as saved:
        time = saved["time"]
        states = saved["state"]
        forecasts = saved["forecast"]
        matlab_command = saved["matlab_command"]
        python_command = saved["python_command"]
        failures = saved["failures"]
    np.testing.assert_allclose(time, np.arange(205, 400.0001, .5),
                               rtol=0, atol=1e-12)
    assert states.shape == (391, 9)
    assert forecasts.shape == (391, 31)
    assert all(array.shape == (391,) for array in (
        matlab_command, python_command, failures,
    ))
    assert all(np.isfinite(array).all() for array in (
        states, forecasts, matlab_command, python_command, failures,
    ))
    assert not failures.any()
    assert np.max(np.abs(matlab_command)) > 1.49e6
    assert np.count_nonzero(np.abs(matlab_command) > 1e3) > 380

    matrices = build_sphere_mpc_matrices(
        MODEL / "sphere.h5", MODEL / "coeff.mat",
    )
    commands = np.empty(len(time))
    for index, (state, forecast) in enumerate(zip(states, forecasts)):
        plan, feasible = solve_sphere_mpc_qp(matrices, state, forecast)
        assert feasible, f"Python QP infeasible at {time[index]:g} s"
        commands[index] = plan[0]
    difference = commands - matlab_command
    assert np.max(np.abs(difference)) < 2.0
    assert np.sqrt(np.mean(difference**2)) < .5

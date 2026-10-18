"""Compare the published WaveStar radiation poles at the source solver step.

This is a one-step solver diagnostic, not a coupled WaveStar trajectory gate.
The amplification factors were measured with Simulink ode8 in R2025b.
"""

from pathlib import Path

import h5py
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import linear_sum_assignment

from wecsim.dynamics.rkdp87 import A as ODE8_A, B as ODE8_B


ROOT = Path(__file__).parent / "fixtures/matlab"
HYDRO = ROOT / "wavestar_fixed_input_step/wavestar.h5"
SOURCE = ROOT / "wavestar_ode8_modes/amplification.npz"


def _active_radiation_poles():
    with h5py.File(HYDRO) as file:
        fit = file["body1/hydro_coeffs/radiation_damping/state_space"]
        order = np.asarray(fit["it"], dtype=int)
        matrix = np.asarray(fit["A/all"])
    # The planar WaveStar linkage excites the float's surge, heave, and pitch.
    return np.concatenate([
        np.linalg.eigvals(matrix[output, motion, :size, :size])
        for output in range(6) for motion in (0, 2, 4)
        if (size := int(order[output, motion])) > 0
    ])


def _single_dop853_step(pole, step):
    solved = solve_ivp(lambda _time, state: pole * state,
                       (0, step), [1 + 0j], method="DOP853",
                       first_step=step, max_step=step,
                       rtol=100, atol=100)
    assert solved.success and solved.t.shape == (2,)
    return solved.y[0, -1]


def test_source_ode8_and_python_step_amplification_on_pinned_modes():
    with np.load(SOURCE, allow_pickle=False) as saved:
        poles = saved["eigenvalues"]
        source_ode8 = saved["ode8_amplification"]
    assert poles.shape == source_ode8.shape == (92,)
    assert np.isfinite(poles).all() and np.isfinite(source_ode8).all()

    # HDF5 is the input; the archived eigenvalues merely pair MATLAB's
    # amplification factor with each pole across BLAS eigenvalue orderings.
    distances = np.abs(_active_radiation_poles()[:, None] - poles[None, :])
    rows, columns = linear_sum_assignment(distances)
    assert np.max(distances[rows, columns]) < 1e-9

    published_step = .05
    z = poles * published_step
    rk4 = 1 + z + z**2 / 2 + z**3 / 6 + z**4 / 24
    assert np.count_nonzero(np.abs(rk4) > 1) == 18
    assert np.max(np.abs(rk4)) > 8.69
    assert np.max(np.abs(source_ode8)) < 1

    fine_z = poles * (published_step / 2)
    fine_rk4 = 1 + fine_z + fine_z**2 / 2 + fine_z**3 / 6 + fine_z**4 / 24
    assert np.max(np.abs(fine_rk4)) < 1

    # SciPy's DOP853 also damps these poles at 50 ms, but it uses a different
    # RK tableau. A single forced step is enough to reject it as an exact
    # replacement for the source's ode8 update.
    dop853 = np.array([_single_dop853_step(pole, published_step)
                       for pole in poles])
    assert np.max(np.abs(dop853)) < 1
    gap = np.max(np.abs(dop853 - source_ode8))
    assert .04 < gap < .043

    # Even exact linear evolution differs from ode8's one-step update, which
    # bounds what a differently integrated Python plant can match exactly.
    exact_gap = np.max(np.abs(np.exp(z) - source_ode8))
    assert .012 < exact_gap < .014


def test_source_ode8_tableau_matches_all_pinned_radiation_modes():
    with np.load(SOURCE, allow_pickle=False) as saved:
        poles = saved["eigenvalues"]
        source_ode8 = saved["ode8_amplification"]
    scaled_poles = poles * .05
    stages = np.zeros((13, len(poles)), dtype=complex)
    for index in range(13):
        stages[index] = scaled_poles * (
            1 + ODE8_A[index, :index] @ stages[:index]
        )
    amplification = 1 + ODE8_B @ stages
    assert np.max(np.abs(amplification - source_ode8)) < 1e-10

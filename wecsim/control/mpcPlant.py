"""Prediction matrices for the published Sphere heave MPC controller."""

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.io import loadmat
from scipy.linalg import expm

from wecsim.hydrodynamics.bodyClass import BodyClass


@dataclass(frozen=True)
class MPCMatrices:
    A: np.ndarray
    Bu: np.ndarray
    Bv: np.ndarray
    C: np.ndarray
    Sx: np.ndarray
    Su: np.ndarray
    Sv: np.ndarray
    Q: np.ndarray
    H: np.ndarray


def build_sphere_mpc_matrices(
    hydro_file: str | Path, coefficient_file: str | Path, *,
    prediction_step: float = .5, prediction_horizon: float = 15,
    rate_penalty: float = 1e-7, rho: float = 1000, g: float = 9.81,
) -> MPCMatrices:
    """Build the source's augmented nine-state heave prediction model.

    The coefficient MAT file supplies the fourth-order radiation fit used by
    the published controller. This function constructs its plant and QP
    matrices; it does not apply a control law to a WEC trajectory.
    """
    count = round(prediction_horizon / prediction_step)
    if (prediction_step <= 0 or count < 1 or
            not np.isclose(count * prediction_step, prediction_horizon,
                           rtol=0, atol=1e-10) or rate_penalty < 0):
        raise ValueError("MPC horizon must contain positive whole prediction steps")
    body = BodyClass(str(hydro_file))
    body.bodyNumber = body.bodyTotal = 1
    body.readH5file()
    if int(np.asarray(body.dof).item()) != 6:
        raise ValueError("Sphere MPC needs a six-DOF hydrodynamic body")
    hydro = body.hydroData["hydro_coeffs"]
    mass = rho * float(np.asarray(body.dispVol).item())
    infinite_added_mass = rho * float(np.asarray(
        hydro["added_mass"]["inf_freq"])[2, 2])
    restoring = rho * g * float(np.asarray(
        hydro["linear_restoring_stiffness"])[2, 2])
    coefficients = loadmat(coefficient_file, simplify_cells=True)["coeff"]
    numerator = np.asarray(coefficients["KradNum"], dtype=float).ravel()
    denominator = np.asarray(coefficients["KradDen"], dtype=float).ravel()
    if (numerator.shape != (4,) or denominator.shape != (5,) or
            not np.isfinite(numerator).all() or
            not np.isfinite(denominator).all() or denominator[0] != 1):
        raise ValueError("MPC coefficient file needs a finite fourth-order radiation fit")
    effective_mass = mass + infinite_added_mass
    if effective_mass <= 0:
        raise ValueError("MPC effective heave mass must be positive")
    A = np.zeros((9, 9))
    A[0, 1] = -restoring / effective_mass
    A[0, 4] = 1 / effective_mass
    A[1:4, 0:3] = np.eye(3)
    A[4, :8] = np.r_[-numerator, -denominator[1:]]
    A[5:8, 4:7] = np.eye(3)
    A[0, 8] = 1 / effective_mass
    Bu = np.r_[np.zeros(8), 1.0]
    Bv = np.r_[1 / effective_mass, np.zeros(8)]
    C = np.zeros((3, 9))
    C[0, 0] = C[1, 1] = C[2, 8] = 1

    # Zero-order hold for the two inputs, matching MATLAB c2d(..., 'zoh').
    augmented = np.zeros((11, 11))
    augmented[:9, :9] = A
    augmented[:9, 9] = Bu
    augmented[:9, 10] = Bv
    discrete = expm(augmented * prediction_step)
    Ad = discrete[:9, :9]
    Bud = discrete[:9, 9]
    Bvd = discrete[:9, 10]
    powers = [np.eye(9)]
    for _ in range(count):
        powers.append(powers[-1] @ Ad)
    Sx = np.vstack([C @ powers[step] for step in range(1, count + 1)])
    Su = np.zeros((3 * count, count + 1))
    Sv = np.zeros_like(Su)
    for step in range(1, count + 1):
        rows = slice(3 * (step - 1), 3 * step)
        for previous in range(step):
            transition = C @ powers[step - previous - 1]
            Su[rows, previous] = transition @ Bud
            Sv[rows, previous] = transition @ Bvd
    q_step = np.array([[0, 0, 1], [0, 0, 0], [1, 0, 0]])
    Q = prediction_step / 2 * np.kron(np.eye(count), q_step)
    H = Su.T @ Q @ Su + rate_penalty * np.eye(count + 1)
    return MPCMatrices(A, Bu, Bv, C, Sx, Su, Sv, Q, H)


def replay_sphere_mpc_plant(
    matrices: MPCMatrices, pto_force: np.ndarray,
    excitation_force: np.ndarray, *, dt: float,
) -> np.ndarray:
    """Replay the controller's internal plant from logged force inputs.

    The result has velocity, position relative to the plant equilibrium,
    and PTO force columns. This is a prediction-model diagnostic, not the
    physical WEC body trajectory or a closed-loop controller.
    """
    force = np.asarray(pto_force, dtype=float)
    excitation = np.asarray(excitation_force, dtype=float)
    if (force.ndim != 1 or excitation.shape != force.shape or len(force) < 2
            or not np.isfinite(force).all() or not np.isfinite(excitation).all()
            or not np.isfinite(dt) or dt <= 0):
        raise ValueError("plant replay needs equal finite force histories and positive dt")
    augmented = np.zeros((11, 11))
    augmented[:9, :9] = matrices.A
    augmented[:9, 9] = matrices.Bu
    augmented[:9, 10] = matrices.Bv
    discrete = expm(augmented * dt)
    Ad = discrete[:9, :9]
    Bud = discrete[:9, 9]
    Bvd = discrete[:9, 10]
    state = np.zeros(9)
    response = np.zeros((len(force), 3))
    for step in range(len(force) - 1):
        rate = (force[step + 1] - force[step]) / dt
        midpoint_excitation = (excitation[step] + excitation[step + 1]) / 2
        state = Ad @ state + Bud * rate + Bvd * midpoint_excitation
        response[step + 1] = matrices.C @ state
    return response

"""Forecasting and force-rate optimization for the published Sphere MPC."""

from dataclasses import dataclass

import numpy as np
from scipy.linalg import expm
from scipy.optimize import Bounds, LinearConstraint, minimize

from wecsim.control.mpcPlant import MPCMatrices


@dataclass(frozen=True)
class MPCResponse:
    time: np.ndarray
    internal_state: np.ndarray
    command_rate: np.ndarray
    pto_force: np.ndarray
    feasible: np.ndarray


def predict_sphere_excitation(
    recent_force: np.ndarray, *, history_steps: int = 200,
    order: int = 4, horizon_steps: int = 30,
) -> np.ndarray:
    """Fit WEC-Sim's autoregression and return current plus future force.

    ``recent_force`` is ordered oldest to newest and includes the current
    sample. The published model uses 205 samples, 201 regression rows,
    four autoregressive coefficients, and 30 predicted samples.
    """
    samples = np.asarray(recent_force, dtype=float)
    required = history_steps + order + 1
    if (samples.ndim != 1 or len(samples) < required
            or not np.isfinite(samples).all()
            or history_steps < order or order < 1 or horizon_steps < 1):
        raise ValueError("MPC forecast needs a finite ordered force history")
    samples = samples[-required:]
    targets = samples[order:][::-1]
    regressors = np.column_stack([
        samples[order - lag - 1:len(samples) - lag - 1][::-1]
        for lag in range(order)
    ])
    coefficients = np.linalg.lstsq(regressors, targets, rcond=None)[0]
    predicted = np.empty(horizon_steps + 1)
    predicted[0] = samples[-1]
    last = list(samples[-order:])
    for step in range(horizon_steps):
        future = float(coefficients @ np.asarray(last[::-1]))
        predicted[step + 1] = future
        last = last[1:] + [future]
    return predicted


def solve_sphere_mpc_qp(
    matrices: MPCMatrices, state: np.ndarray,
    predicted_excitation: np.ndarray, *,
    max_speed: float = 3.0, max_position: float = 4.0,
    max_force: float = 2e6, max_force_rate: float = 1.5e6,
) -> tuple[np.ndarray, bool]:
    """Solve the published force-rate QP from one sampled plant state.

    Returns the full force-rate plan and whether the constrained solve
    succeeded. As in ``mpcFcn.m``, an infeasible solve supplies zero rate.
    The first plan entry is the command applied to the PTO force integrator.
    """
    x = np.asarray(state, dtype=float)
    forecast = np.asarray(predicted_excitation, dtype=float)
    count = matrices.Su.shape[1]
    limits = np.asarray([max_speed, max_position, max_force,
                         max_force_rate], dtype=float)
    if (x.shape != (matrices.A.shape[0],)
            or forecast.shape != (count,)
            or not np.isfinite(x).all() or not np.isfinite(forecast).all()
            or not np.isfinite(limits).all() or np.any(limits <= 0)):
        raise ValueError("MPC QP needs finite state, forecast, and positive limits")
    baseline = matrices.Sx @ x + matrices.Sv @ forecast
    gradient = matrices.Su.T @ matrices.Q @ baseline
    output_limit = np.tile(limits[:3], count - 1)
    output_scale = output_limit
    # Numerical scaling changes neither objective minimizer nor constraints.
    force_scale = 1e6
    objective_scale = 1e6
    hessian = matrices.H * force_scale**2 / objective_scale
    linear = gradient * force_scale / objective_scale
    output_matrix = matrices.Su * force_scale / output_scale[:, None]
    lower = (-output_limit - baseline) / output_scale
    upper = (output_limit - baseline) / output_scale
    solution = minimize(
        lambda rate: .5 * rate @ hessian @ rate + linear @ rate,
        np.zeros(count),
        jac=lambda rate: hessian @ rate + linear,
        method="SLSQP",
        bounds=Bounds(-max_force_rate / force_scale * np.ones(count),
                      max_force_rate / force_scale * np.ones(count)),
        constraints=[LinearConstraint(output_matrix, lower, upper)],
        options={"ftol": 1e-13, "maxiter": 600},
    )
    if not solution.success:
        return np.zeros(count), False
    return solution.x * force_scale, True


def integrate_sphere_mpc_force(
    command_rate: np.ndarray, *, dt: float, end_time: float,
    control_step: float = .5, command_delay: float = .5,
) -> np.ndarray:
    """Integrate held MPC force-rate commands through the source delay."""
    commands = np.asarray(command_rate, dtype=float)
    if (commands.ndim != 1 or not np.isfinite(commands).all()
            or not np.isfinite([dt, end_time, control_step, command_delay]).all()
            or dt <= 0 or end_time < 0 or control_step <= 0
            or command_delay < 0):
        raise ValueError("MPC force integration needs finite rates and times")
    steps = round(end_time / dt)
    stride = round(control_step / dt)
    delay = round(command_delay / dt)
    if (stride < 1 or not np.isclose(steps * dt, end_time, atol=1e-10)
            or not np.isclose(stride * dt, control_step, atol=1e-10)
            or not np.isclose(delay * dt, command_delay, atol=1e-10)
            or len(commands) < steps // stride + 1):
        raise ValueError("MPC timing must align to dt and cover the run")
    sampled_index = (np.arange(steps) - delay) // stride
    applied_rate = np.zeros(steps)
    valid = sampled_index >= 0
    applied_rate[valid] = commands[sampled_index[valid]]
    return np.r_[0.0, np.cumsum(applied_rate * dt)]


def simulate_sphere_mpc(
    matrices: MPCMatrices, excitation_force: np.ndarray, *,
    dt: float = .01, control_step: float = .5,
    command_delay: float = .5, start_time: float = 205,
    history_steps: int = 200, order: int = 4,
    max_speed: float = 3.0, max_position: float = 4.0,
    max_force: float = 2e6, max_force_rate: float = 1.5e6,
) -> MPCResponse:
    """Run the published sampled MPC against an independent force history.

    The prediction plant and optimizer use wave excitation, not physical
    body displacement. The published Simulink feedback advances its start
    counter every other 0.5 s controller sample, making 205 s its first
    active command; callers may choose a different ``start_time``.
    """
    excitation = np.asarray(excitation_force, dtype=float)
    if (excitation.ndim != 1 or len(excitation) < 2
            or not np.isfinite(excitation).all()
            or not np.isfinite([dt, control_step, command_delay, start_time]).all()
            or dt <= 0 or control_step <= 0 or command_delay < 0
            or start_time < 0 or history_steps < order or order < 1):
        raise ValueError("Sphere MPC needs finite excitation and aligned positive timing")
    stride = round(control_step / dt)
    delay = round(command_delay / dt)
    start = round(start_time / dt)
    steps = len(excitation) - 1
    if (stride < 1 or not np.isclose(stride * dt, control_step, atol=1e-10)
            or not np.isclose(delay * dt, command_delay, atol=1e-10)
            or not np.isclose(start * dt, start_time, atol=1e-10)
            or start % stride or start < history_steps * stride):
        raise ValueError("MPC start and timing must align with a full force history")
    count = steps // stride + 1
    commands = np.zeros(count)
    feasible = np.ones(count, dtype=bool)
    state = np.zeros((len(excitation), matrices.A.shape[0]))
    augmented = np.zeros((11, 11))
    augmented[:9, :9] = matrices.A
    augmented[:9, 9] = matrices.Bu
    augmented[:9, 10] = matrices.Bv
    discrete = expm(augmented * dt)
    Ad = discrete[:9, :9]
    Bud = discrete[:9, 9]
    Bvd = discrete[:9, 10]
    for step in range(len(excitation)):
        if step % stride == 0 and step >= start:
            command_index = step // stride
            recent = excitation[
                step - (history_steps + order) * stride:step + 1:stride
            ]
            forecast = predict_sphere_excitation(
                recent, history_steps=history_steps, order=order,
                horizon_steps=matrices.Su.shape[1] - 1,
            )
            plan, feasible[command_index] = solve_sphere_mpc_qp(
                matrices, state[step], forecast,
                max_speed=max_speed, max_position=max_position,
                max_force=max_force, max_force_rate=max_force_rate,
            )
            commands[command_index] = plan[0]
        if step == steps:
            break
        held_index = (step - delay) // stride
        held_rate = commands[held_index] if held_index >= 0 else 0.0
        midpoint_excitation = (excitation[step] + excitation[step + 1]) / 2
        state[step + 1] = (Ad @ state[step] + Bud * held_rate
                           + Bvd * midpoint_excitation)
    time = np.arange(len(excitation)) * dt
    return MPCResponse(time, state, commands, state[:, 8].copy(), feasible)

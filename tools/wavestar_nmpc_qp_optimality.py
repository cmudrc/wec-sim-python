"""Check the pinned WaveStar RTI QP at selected sensitive source samples.

The source state, forecast, and preceding plan are fixed inputs to this
component diagnostic. No result here constitutes autonomous motion parity.
"""

import argparse
import json
import platform
from hashlib import sha256
from pathlib import Path

import numpy as np
import scipy
from scipy.io import loadmat

import wecsim.control.wavestarNmpc as nmpc


ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "tests/fixtures/matlab/wavestar_nmpc_reference_convergence_225s"
MATRICES = ROOT / "tests/fixtures/matlab/wavestar_nmpc_full_step/direct_controller_inputs.npz"
HASHES = {
    "controller_0p25ms_compatible_amd7763.npz":
        "06d55c1972b6e9130fbeeea52d29cc6691fddfb7408c687c950ce9522c2498e3",
    "checked_qp_replay_compatible_amd7763.mat":
        "cbc5a2acfe651f3570db0405d040ee3d40791b07f30ac4b3d282ae289ec5af7d",
    "direct_controller_inputs.npz":
        "3ff8f4ba8654ef66bf9058a9fe717cfa6c60868300d19e64a089fcb06781dbf0",
}


def _checked(path):
    digest = sha256(path.read_bytes()).hexdigest()
    if digest != HASHES[path.name]:
        raise ValueError(f"Pinned QP input changed: {path}: {digest}")


def _qp_metrics(hessian, linear, bounds, delta):
    gradient = hessian @ delta + linear
    lower, upper = bounds[:, 0], bounds[:, 1]
    at_lower = delta <= lower + 1e-9
    at_upper = delta >= upper - 1e-9
    free = ~(at_lower | at_upper)
    residual = np.zeros_like(delta)
    residual[at_lower] = np.minimum(gradient[at_lower], 0.)
    residual[at_upper] = np.maximum(gradient[at_upper], 0.)
    residual[free] = gradient[free]
    return {
        "qp_objective": float(.5 * delta @ hessian @ delta + linear @ delta),
        "max_projected_kkt_residual": float(np.max(np.abs(residual))),
        "max_bound_violation_Nm": float(max(
            0., np.max(lower - delta), np.max(delta - upper),
        )),
        "free_variables_at_1e-9_Nm_tolerance": int(np.count_nonzero(free)),
    }


def _sample(time_s, source_prior, states, forecasts, plans, matrices):
    index = round(time_s / .05)
    active = index - 300
    controller = nmpc.WaveStarNmpcController()
    controller._transition = matrices["A"].copy()
    controller._input = matrices["Bu"][:, 0].copy()
    controller._wave_input = matrices["Bw"][:, 0].copy()
    if source_prior and active > 0:
        controller._nominal = plans[active - 1].copy()
    shifted = np.r_[controller._nominal[1:], controller._nominal[-1]]
    captured = {}
    original = nmpc.minimize

    def capture(objective, x0, *, jac, bounds, method, options):
        if method == "L-BFGS-B":
            zero = np.zeros(40)
            linear = jac(zero)
            basis = np.eye(40)
            hessian = np.column_stack([
                jac(basis[i]) - linear for i in range(40)
            ])
            captured.update(hessian=hessian, linear=linear,
                            bounds=np.asarray(bounds))
        return original(objective, x0, jac=jac, bounds=bounds,
                        method=method, options=options)

    nmpc.minimize = capture
    try:
        controller.step(states[index], forecasts[index], time_s)
    finally:
        nmpc.minimize = original

    hessian = captured["hessian"]
    linear = captured["linear"]
    bounds = captured["bounds"]
    python_delta = controller._nominal - shifted
    source_delta = plans[active] - shifted
    return {
        "time_s": time_s,
        "prior_plan": "source" if source_prior else "zero",
        "hessian_min_eigenvalue": float(np.linalg.eigvalsh(hessian).min()),
        "hessian_condition_number": float(np.linalg.cond(hessian)),
        "python_first_request_Nm": float(controller._nominal[0]),
        "matlab_first_request_Nm": float(plans[active, 0]),
        "python_solution": _qp_metrics(hessian, linear, bounds, python_delta),
        "matlab_plan_in_same_qp": _qp_metrics(
            hessian, linear, bounds, source_delta,
        ),
    }


def audit():
    controller_file = REF / "controller_0p25ms_compatible_amd7763.npz"
    replay_file = REF / "checked_qp_replay_compatible_amd7763.mat"
    for path in (controller_file, replay_file, MATRICES):
        _checked(path)
    with np.load(controller_file, allow_pickle=False) as saved:
        time = saved["time"].copy()
        states = saved["estimated_state"].copy()
        forecasts = saved["forecast"].copy()
    with np.load(MATRICES, allow_pickle=False) as saved:
        matrices = {key: saved[key].copy() for key in ("A", "Bu", "Bw")}
    replay = loadmat(replay_file)
    plans = replay["plans"]
    np.testing.assert_allclose(time, np.arange(4501) * .05,
                               rtol=0, atol=1e-9)
    np.testing.assert_array_equal(replay["time"].ravel(), time[300:])
    if (states.shape != (4501, 5) or forecasts.shape != (4501, 40) or
            plans.shape != (4201, 40)):
        raise ValueError("Pinned QP arrays changed shape")
    cases = [
        _sample(t, True, states, forecasts, plans, matrices)
        for t in (17.5, 30., 97.25, 202.25)
    ]
    cases.append(_sample(30., False, states, forecasts, plans, matrices))
    return {
        "purpose": "fixed-input RTI QP optimality and plan-policy diagnostic",
        "pinned_sha256": HASHES,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "active_bound_tolerance_Nm": 1e-9,
        "samples": cases,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = audit()
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()

"""Check whether coarser source logging changes ellipsoid ode45 dynamics.

Both copied R2025b inputs use the same 0.00125 s solver, nonlinear-force,
radiation-convolution, and maximum-step clocks. Only `simu.dtOut` differs.
This diagnostic does not replace the all-clock source convergence gate.
"""

import os
from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import trapezoid


RECORDS = os.environ.get("WEC_SIM_ELLIPSOID_AB_DIR")
CASE = os.environ.get("WEC_SIM_ELLIPSOID_AB_CASE")
pytestmark = pytest.mark.skipif(
    not (RECORDS and CASE),
    reason="fresh pinned R2025b output-grid A/B output absent",
)

# This A/B is stricter than the derived source final-step motion/PTO screen.
AB_LIMITS = {
    1: (1e-5, "heave position"),
    2: (2e-5, "heave speed"),
    3: (1e-5, "PTO stroke"),
    4: (2e-5, "PTO speed"),
    5: (30, "PTO force"),
    6: (30, "absorbed power"),
}
FULL_TO_SAMPLED_WORK_LIMIT = 1e-5
AB_WORK_LIMIT = 1e-4


def _read(label, output_step):
    folder = Path(RECORDS) / CASE
    history = np.loadtxt(folder / f"{label}_history.csv",
                         delimiter=",", ndmin=2)
    settings = np.loadtxt(folder / f"{label}_settings.csv",
                          delimiter=",", ndmin=2)
    assert history.shape == (round(150 / output_step) + 1, 11)
    assert np.isfinite(history).all()
    np.testing.assert_allclose(settings,
                               [[.00125, output_step, .00125, .00125, .00125]],
                               rtol=0, atol=1e-12)
    np.testing.assert_allclose(history[:, 0],
                               np.arange(len(history)) * output_step,
                               rtol=0, atol=1e-9)
    return history


def _maximum(a, b):
    assert a.shape == b.shape
    return float(np.max(np.abs(a - b)))


def test_pinned_output_only_ab_preserves_motion_pto_and_work():
    assert CASE in ("Regular", "RegularCIC")
    full = _read("all125", .00125)
    output = _read("out005", .005)
    common = full[::4]
    assert _maximum(common[:, 0], output[:, 0]) < 1e-9
    assert _maximum(common[:, 7], output[:, 7]) < 1e-10
    for column, (limit, label) in AB_LIMITS.items():
        error = _maximum(common[:, column], output[:, column])
        print(f"output-only A/B {label}: {error:.8g} (gate {limit:.8g})")
        assert error < limit
    for column, label in ((8, "restoring"), (9, "excitation"),
                          (10, "acceleration")):
        print(f"output-only A/B {label} diagnostic:",
              _maximum(common[:, column], output[:, column]))

    full_work = float(trapezoid(full[:, 6], full[:, 0]))
    sampled_work = float(trapezoid(common[:, 6], common[:, 0]))
    output_work = float(trapezoid(output[:, 6], output[:, 0]))
    assert min(full_work, sampled_work, output_work) > 1e3
    sampling_error = abs(sampled_work / full_work - 1)
    ab_error = abs(output_work / sampled_work - 1)
    print("signed PTO work (J):", full_work, sampled_work, output_work)
    print("full/sampled signed-work change:", sampling_error,
          "(gate", FULL_TO_SAMPLED_WORK_LIMIT, ")")
    print("output-only A/B signed-work change:", ab_error,
          "(gate", AB_WORK_LIMIT, ")")
    assert sampling_error < FULL_TO_SAMPLED_WORK_LIMIT
    assert ab_error < AB_WORK_LIMIT

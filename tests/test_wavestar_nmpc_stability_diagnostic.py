"""Preserve the pinned evidence for WaveStar's open pointwise gate."""

import numpy as np
import pytest

from tools.wavestar_nmpc_autonomous_step_sensitivity import (
    LIMITS, _compare, _load_source,
)
from tools.wavestar_nmpc_qp_optimality import audit as qp_audit


class PublishedClosedLoopStepGap(AssertionError):
    """The checked source pair exceeds an unchanged pointwise limit."""


@pytest.mark.xfail(
    strict=True,
    reason="published RTI closed loop changes beyond the original pointwise gate",
    raises=PublishedClosedLoopStepGap,
)
def test_pinned_matlab_source_step_stability_under_original_limits():
    half = _load_source("source_0p125ms_macos_arm_m1.npz")
    quarter = _load_source("source_0p0625ms_macos_arm_m1.npz")
    np.testing.assert_array_equal(half["components"], quarter["components"])
    comparison = _compare(half, quarter)
    failed = {
        label: comparison["pointwise"][label]["max_abs"]
        for label, limit in LIMITS.items()
        if comparison["pointwise"][label]["max_abs"] >= limit
    }
    if failed:
        raise PublishedClosedLoopStepGap(str(failed))


def test_pinned_qp_is_well_conditioned_but_plan_policy_changes_objective():
    report = qp_audit()
    source_prior = [sample for sample in report["samples"]
                    if sample["prior_plan"] == "source"]
    assert len(source_prior) == 4
    for sample in source_prior:
        assert sample["hessian_min_eigenvalue"] > .8
        assert sample["hessian_condition_number"] < 2
        assert sample["python_solution"]["max_bound_violation_Nm"] == 0
        assert (abs(sample["python_first_request_Nm"] -
                    sample["matlab_first_request_Nm"]) < 2e-7)
    reset = report["samples"][-1]
    assert reset["time_s"] == 30. and reset["prior_plan"] == "zero"
    assert (abs(reset["python_first_request_Nm"] -
                reset["matlab_first_request_Nm"]) > 1.)

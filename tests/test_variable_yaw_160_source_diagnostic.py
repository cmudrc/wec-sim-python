"""Keep the local 160 s MATLAB selected-heading step limit inspectable."""

import json
from pathlib import Path

import numpy as np


FIXTURES = Path(__file__).parent / "fixtures/matlab"


def test_matlab_quarter_step_loses_stability_between_140_and_160_seconds():
    with np.load(FIXTURES / "oswec_yaw_160_source.npz") as saved:
        fine = saved["fine"]
        quarter = saved["quarter"]
        fine_wave = saved["fine_wave"]
        quarter_wave = saved["quarter_wave"]
    reported = json.loads((FIXTURES / "oswec_yaw_160_metrics.json").read_text())
    gates = json.loads((FIXTURES / "oswec_yaw_160_gate_spec.json").read_text())
    limits = tuple(gates["source_self_change_limits"][key] for key in
                   ("yaw_rad", "yaw_speed_rad_s", "pto_torque_Nm",
                    "relative_pto_work"))
    assert gates["source_self_change_limits"] == gates["independent_python_source_limits"]
    assert fine.shape == (50_001, 7)
    assert quarter.shape == (64_001, 7)
    assert fine_wave.shape == (50_001, 2)
    assert quarter_wave.shape == (64_001, 2)
    assert all(np.isfinite(array).all() for array in
               (fine, quarter, fine_wave, quarter_wave))
    np.testing.assert_allclose(quarter[::2, 0], fine[:32_001, 0], rtol=0, atol=1e-10)
    np.testing.assert_allclose(quarter_wave[::2], fine_wave[:32_001],
                               rtol=0, atol=1e-12)
    np.testing.assert_allclose(quarter[:, 5], -120_000 * quarter[:, 2],
                               rtol=0, atol=1e-8)
    np.testing.assert_allclose(quarter[:, 6], quarter[:, 5] * quarter[:, 2],
                               rtol=0, atol=1e-8)

    # These limits were fixed for the earlier 120 and 140 s source audits.
    # A green diagnostic here records the 160 s failure; it is not a parity gate.
    for duration, stable in ((120, True), (140, True), (160, False)):
        q = quarter[:round(duration / .0025) + 1]
        f = fine[:round(duration / .005) + 1]
        q_work = -np.trapezoid(q[:, 6], q[:, 0])
        f_work = -np.trapezoid(f[:, 6], f[:, 0])
        measured = (
            np.max(np.abs(f[:, 1] - q[::2, 1])),
            np.max(np.abs(f[:, 2] - q[::2, 2])),
            np.max(np.abs(f[:, 5] - q[::2, 5])),
            abs(q_work - f_work) / q_work,
        )
        np.testing.assert_allclose(measured, reported[str(duration)]["source_step"],
                                   rtol=0, atol=1e-12)
        assert all(value < limit for value, limit in zip(measured, limits)) is stable
        assert np.count_nonzero(f[:, 3] != q[::2, 3]) == (
            reported[str(duration)]["source_heading_mismatch_samples"]
        )

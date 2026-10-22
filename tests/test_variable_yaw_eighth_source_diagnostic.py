"""Preserve the nonmonotone 160 s MATLAB selected-heading step result."""

import hashlib
import gzip
import json
from pathlib import Path

import numpy as np


FIXTURES = Path(__file__).parent / "fixtures/matlab"
CHANNELS = ("yaw_rad", "yaw_speed_rad_s", "pto_torque_Nm", "relative_pto_work")


def _errors(reference, actual):
    reference_work = -np.trapezoid(reference[:, 6], reference[:, 0])
    actual_work = -np.trapezoid(actual[:, 6], actual[:, 0])
    assert reference_work > 0 and actual_work > 0
    return {
        "yaw_rad": np.max(np.abs(actual[:, 1] - reference[:, 1])),
        "yaw_speed_rad_s": np.max(np.abs(actual[:, 2] - reference[:, 2])),
        "pto_torque_Nm": np.max(np.abs(actual[:, 5] - reference[:, 5])),
        "relative_pto_work": abs(actual_work - reference_work) / reference_work,
    }


def test_matlab_eighth_step_branch_and_unchanged_bounded_gates():
    archive = FIXTURES / "oswec_yaw_160_eighth_source.npz"
    evidence = json.loads((FIXTURES / "oswec_yaw_160_eighth_evidence.json").read_text())
    gate_file = FIXTURES / "oswec_yaw_160_gate_spec.json"
    gates = json.loads(gate_file.read_text())
    assert hashlib.sha256(archive.read_bytes()).hexdigest() == evidence["eighth_fixture_sha256"]
    assert hashlib.sha256(gate_file.read_bytes()).hexdigest() == evidence["gate_spec_sha256"]
    runner = FIXTURES / "oswec_yaw_160_eighth_runner.m"
    assert hashlib.sha256(runner.read_bytes()).hexdigest() == evidence["source_runner_sha256"]
    source_input = FIXTURES / "oswec_yaw_160_eighth_input.m.gz"
    assert hashlib.sha256(gzip.decompress(source_input.read_bytes())).hexdigest() == evidence[
        "source_input_sha256"
    ]
    assert gates["published_case_promotion"] is evidence["published_case_promotion"] is False
    with np.load(archive, allow_pickle=False) as saved:
        eighth = saved["eighth"]
        eighth_wave = saved["eighth_wave"]
        components = saved["components"]
        quarter_components = saved["quarter_components"]
    with np.load(FIXTURES / "oswec_yaw_160_source.npz", allow_pickle=False) as saved:
        half = saved["fine"][:32_001]
        quarter = saved["quarter"]
        quarter_wave = saved["quarter_wave"]
    assert eighth.shape == (128_001, 7)
    assert eighth_wave.shape == (128_001, 2)
    assert components.shape == (500, 4)
    assert quarter_components.shape == (500, 4)
    assert all(np.isfinite(a).all() for a in
               (eighth, eighth_wave, components, quarter_components))
    np.testing.assert_array_equal(components, quarter_components)
    np.testing.assert_allclose(eighth[::2, 0], quarter[:, 0], rtol=0, atol=1e-10)
    np.testing.assert_allclose(eighth_wave[::2], quarter_wave, rtol=0, atol=1e-12)
    np.testing.assert_allclose(eighth[:, 5], -120_000 * eighth[:, 2],
                               rtol=0, atol=1e-8)
    np.testing.assert_allclose(eighth[:, 6], eighth[:, 5] * eighth[:, 2],
                               rtol=0, atol=1e-8)
    assert not np.any((eighth[:, 4] <= -12) | (eighth[:, 4] >= 12))

    limits = gates["source_self_change_limits"]
    for duration, row in zip((120, 140, 160), evidence["paired_metrics"]["windows"],
                             strict=True):
        q = quarter[:round(duration / .0025) + 1]
        e = eighth[:round(duration / .00125) + 1]
        measured = _errors(q, e[::2])
        # Work uses every eighth-step sample, including those between quarter steps.
        measured["relative_pto_work"] = abs(
            -np.trapezoid(e[:, 6], e[:, 0])
            + np.trapezoid(q[:, 6], q[:, 0])
        ) / -np.trapezoid(q[:, 6], q[:, 0])
        for channel in CHANNELS:
            np.testing.assert_allclose(measured[channel], row["source_self"][channel],
                                       rtol=0, atol=1e-12)
            assert bool(measured[channel] < limits[channel]) is (duration < 160)
        assert row["source_self_pass"] is (duration < 160)
        assert np.count_nonzero(e[::2, 4] != q[:, 4]) == row[
            "source_heading_mismatch_samples"
        ]

    # The quarter step briefly enters the -3 degree bank at the local peak;
    # both neighboring source steps remain below the same selection threshold.
    threshold = np.deg2rad(12.875)
    for trace, crosses in ((half, False), (quarter, True), (eighth, False)):
        local = trace[(trace[:, 0] >= 144) & (trace[:, 0] <= 146)]
        assert bool(np.max(local[:, 1]) > threshold) is crosses
    quarter_events = np.flatnonzero(np.diff(quarter[:, 4])) + 1
    eighth_events = np.flatnonzero(np.diff(eighth[:, 4])) + 1
    assert np.any((quarter[quarter_events, 0] == 144.8475)
                  & (quarter[quarter_events, 4] == -3))
    assert np.any((eighth[eighth_events, 0] == 146.45)
                  & (eighth[eighth_events, 4] == -3))
    assert len(quarter_events) == 257
    assert len(eighth_events) == 256

    # Half and eighth steps are close, but the intervening quarter step is not:
    # this does not justify a convergence or published-case parity claim.
    half_to_eighth = _errors(half, eighth[::4])
    half_to_eighth["relative_pto_work"] = abs(
        -np.trapezoid(eighth[:, 6], eighth[:, 0])
        + np.trapezoid(half[:, 6], half[:, 0])
    ) / -np.trapezoid(half[:, 6], half[:, 0])
    assert all(half_to_eighth[channel] < limits[channel] for channel in CHANNELS)
    assert evidence["paired_metrics"]["all_bounded_gates_pass"] is False

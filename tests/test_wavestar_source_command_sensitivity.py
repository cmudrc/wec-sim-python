"""Measure the published WaveStar controller's source-side feedback sensitivity."""

from pathlib import Path

import numpy as np


SOURCE = (Path(__file__).parent / "fixtures/matlab/wavestar_nmpc_full_step"
          / "source_command_pulse.npz")
SOURCE_FIRST_SPLIT = (Path(__file__).parent
                      / "fixtures/matlab/wavestar_nmpc_full_step"
                      / "source_command_pulse_2p7e.npz")


def test_single_small_source_command_pulse_amplifies_in_closed_loop():
    # Both records are fresh local R2025b 225 s runs with the same saved HDF5
    # and incident sea. The second changes one active NMPC request by 1e-4 Nm.
    with np.load(SOURCE, allow_pickle=False) as source:
        time = source["time"]
        baseline_pitch = source["baseline_pitch"]
        baseline_command = source["baseline_command"]
        pulse_pitch = source["pulse_pitch"]
        pulse_command = source["pulse_command"]
    assert all(values.shape == (4501,) for values in (
        time, baseline_pitch, baseline_command, pulse_pitch, pulse_command,
    ))
    assert all(np.isfinite(values).all() for values in (
        time, baseline_pitch, baseline_command, pulse_pitch, pulse_command,
    ))
    np.testing.assert_allclose(time, np.arange(4501) * .05, rtol=0, atol=1e-10)
    pulse_index = 331  # 16.55 s, the 32nd active controller call.
    np.testing.assert_array_equal(
        pulse_command[:pulse_index], baseline_command[:pulse_index],
    )
    np.testing.assert_array_equal(
        pulse_pitch[:pulse_index + 1], baseline_pitch[:pulse_index + 1],
    )
    np.testing.assert_allclose(
        pulse_command[pulse_index] - baseline_command[pulse_index],
        1e-4, rtol=0, atol=1e-12,
    )
    # This describes source sensitivity; it is not a Python parity gate.
    assert np.max(np.abs(pulse_pitch - baseline_pitch)) > .01


def test_first_split_sized_source_pulse_changes_published_response():
    """Measure MATLAB's response to one pulse near the Python first split."""
    with np.load(SOURCE_FIRST_SPLIT, allow_pickle=False) as source:
        values = {name: source[name] for name in source.files}
    assert all(value.shape == (4501,) for value in values.values())
    assert all(np.isfinite(value).all() for value in values.values())
    time = values["time"]
    np.testing.assert_allclose(time, np.arange(4501) * .05,
                               rtol=0, atol=1e-10)
    first = 331  # 32nd active NMPC request, at 16.55 s.
    np.testing.assert_array_equal(
        values["pulse_command"][:first],
        values["baseline_command"][:first],
    )
    np.testing.assert_array_equal(
        values["pulse_pitch"][:first + 1],
        values["baseline_pitch"][:first + 1],
    )
    np.testing.assert_allclose(
        values["pulse_command"][first]
        - values["baseline_command"][first],
        .00272438, rtol=0, atol=1e-10,
    )
    pitch_change = np.max(np.abs(
        values["pulse_pitch"] - values["baseline_pitch"]
    ))
    peak_change = abs(
        np.max(np.abs(values["pulse_pitch"]))
        - np.max(np.abs(values["baseline_pitch"]))
    ) / np.max(np.abs(values["baseline_pitch"]))
    speed_change = np.max(np.abs(
        values["pulse_speed"] - values["baseline_speed"]
    ))
    baseline_work = np.trapezoid(
        values["baseline_pto_force"] * values["baseline_pto_speed"], time,
    )
    pulse_work = np.trapezoid(
        values["pulse_pto_force"] * values["pulse_pto_speed"], time,
    )
    # These describe source self-change, not independent Python parity.
    assert .009 < pitch_change < .011
    assert .007 < peak_change < .009
    assert .05 < speed_change < .06
    assert abs(pulse_work - baseline_work) / abs(baseline_work) < .001

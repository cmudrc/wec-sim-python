"""Measure the published WaveStar controller's source-side feedback sensitivity."""

from pathlib import Path

import numpy as np


SOURCE = (Path(__file__).parent / "fixtures/matlab/wavestar_nmpc_full_step"
          / "source_command_pulse.npz")


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

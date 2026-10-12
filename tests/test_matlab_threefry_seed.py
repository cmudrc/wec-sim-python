"""Check the pinned MATLAB Threefry stream independently of WEC motion."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim import PMWave
from wecsim.environment.irregularWave import _random_phases


@pytest.mark.parametrize("substream", (1, 2, 3))
def test_matlab_threefry_stream(substream):
    directory = os.environ.get("WEC_SIM_MATLAB_THREEFRY_DIR")
    if not directory:
        pytest.skip("live MATLAB Threefry reference not provided")
    reference = np.loadtxt(
        Path(directory) / f"seed1_substream{substream}.csv", delimiter=",")
    generated = _random_phases((1, len(reference)), substream, "matlab")[0]
    np.testing.assert_allclose(generated, 2 * np.pi * reference,
                               rtol=0, atol=1e-14)


def test_numpy_seed_remains_default():
    expected = 2 * np.pi * np.random.default_rng(5).random((3, 2))
    np.testing.assert_array_equal(_random_phases((3, 2), 5, "numpy"), expected)
    assert PMWave(2, 5, seed=5).as_case().get("phase_generator") is None
    assert PMWave(2, 5, seed=5, phase_generator="matlab").as_case()[
        "phase_generator"] == "matlab"


@pytest.mark.parametrize("seed", (None, 0, -1, True, 1.5))
def test_matlab_substream_needs_positive_integer(seed):
    with pytest.raises(ValueError, match="positive integer substream seed"):
        _random_phases((3, 1), seed, "matlab")


def test_unknown_phase_generator_is_rejected():
    with pytest.raises(ValueError, match="phase_generator"):
        PMWave(2, 5, seed=1, phase_generator="unknown").as_case()

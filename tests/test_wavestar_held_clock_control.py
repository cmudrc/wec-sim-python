"""Control: explicit 50 ms controller scheduling must preserve the source run."""

import os
from pathlib import Path

import numpy as np
import pytest


ROOT = Path(__file__).parent
SOURCE = os.environ.get("WEC_SIM_WAVESTAR_CLOCK_CONTROL_DIR")
BASELINE = ROOT / "fixtures/matlab/wavestar_nmpc_full_step/baseline.npz"


def _read(name):
    values = np.loadtxt(Path(SOURCE) / f"{name}.csv", delimiter=",", ndmin=2)
    assert np.isfinite(values).all(), name
    return values


def _gate(name, actual, expected, limit):
    assert actual.shape == expected.shape, name
    error = float(np.max(np.abs(actual - expected)))
    assert error < limit, f"{name}: {error:.9g} >= {limit:.9g}"


@pytest.mark.skipif(not SOURCE, reason="fresh 50 ms MATLAB control run absent")
def test_explicit_controller_clock_preserves_published_50ms_source():
    """Otherwise the 25 ms comparison changes more than the plant step."""
    with np.load(BASELINE, allow_pickle=False) as saved:
        body = saved["body"]
        pto = saved["pto"]
        command = np.column_stack((saved["command_time"], saved["command"]))
        wave = saved["wave"]
        components = saved["components"]
    _gate("incident wave", _read("wave"), wave, 1e-12)
    _gate("wave components", _read("components"), components, 1e-12)
    _gate("body position and velocity", _read("float"), body[:, :13], 1e-8)
    _gate("controller command", _read("command"), command, 1e-6)
    source_axial = _read("pto_axial")
    published_axial = pto[:, [0, 33, 9]]
    _gate("axial PTO force and speed", source_axial, published_axial, 1e-6)
    source_work = np.trapezoid(source_axial[:, 1] * source_axial[:, 2],
                               x=source_axial[:, 0])
    published_work = np.trapezoid(published_axial[:, 1] * published_axial[:, 2],
                                  x=published_axial[:, 0])
    work_change = abs(source_work - published_work) / max(abs(published_work), 1e-12)
    assert work_change < 1e-6, f"relative PTO work change: {work_change:.9g}"

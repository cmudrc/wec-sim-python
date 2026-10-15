"""Measure full-duration source motion when only WaveStar's ode8 step halves.

The earlier 225 s audit stopped before exporting body and PTO histories when
it found that the controller also ran at 25 ms. This diagnostic preserves
that observation and measures the resulting source-to-source motion change;
it does not treat the run as a plant-only parity reference.
"""

import json
import os
from pathlib import Path

import h5py
import numpy as np
import pytest
from scipy.io import loadmat


ROOT = Path(__file__).parent
SOURCE = os.environ.get("WEC_SIM_WAVESTAR_FULL_INTERNAL_STEP_DIR")
HYDRO = os.environ.get("WEC_SIM_WAVESTAR_H5")
PREFIX = "WECCCOMP_NMPC_INTERNAL_STEP_FULL_"
BASELINE = ROOT / "fixtures/matlab/wavestar_nmpc_full_step/baseline.npz"
PRECONTROL = ROOT / "fixtures/matlab/wavestar_fixed_input_step/source.npz"


def _fresh(label):
    values = np.loadtxt(Path(SOURCE) / f"{PREFIX}{label}.csv",
                        delimiter=",", ndmin=2)
    assert np.isfinite(values).all(), label
    return values


def _common_grid(values):
    time = values[:, 0]
    if len(time) == 4501:
        np.testing.assert_allclose(time, np.arange(4501) * .05,
                                   rtol=0, atol=1e-9)
        return values
    assert len(time) == 9001
    np.testing.assert_allclose(time, np.arange(9001) * .025,
                               rtol=0, atol=1e-9)
    return values[::2]


def _difference(a, b):
    delta = np.asarray(a) - np.asarray(b)
    return {"maximum": float(np.max(np.abs(delta))),
            "rms": float(np.sqrt(np.mean(delta**2)))}


@pytest.mark.skipif(not (SOURCE and HYDRO),
                    reason="fresh full-duration MATLAB source run absent")
def test_full_source_self_change_with_controller_clock_observed():
    with np.load(BASELINE) as saved:
        baseline = {key: saved[key] for key in saved.files}
    with np.load(PRECONTROL) as saved:
        short = {key: saved[key] for key in saved.files}
    body_raw = _fresh("WECCCOMP_Nonlinear_Model_Predictive_body1")
    pto_raw = _fresh("WECCCOMP_Nonlinear_Model_Predictive_pto1")
    body = _common_grid(body_raw)
    pto = _common_grid(pto_raw)
    wave = _fresh("wave")
    components = _fresh("components")
    steps = _fresh("steps")
    controller = loadmat(Path(SOURCE) / f"{PREFIX}controller.mat",
                         simplify_cells=True)
    command_time = np.asarray(controller["cmd_ptoM"]["time"]).ravel()
    command = np.asarray(controller["cmd_ptoM"]["signals"]["values"]).ravel()
    power = np.asarray(controller["Output_power"]["signals"]["values"]).ravel()
    energy = np.asarray(controller["Output_energy"]["signals"]["values"]).ravel()

    assert body.shape == (4501, 25)
    assert pto.shape == (4501, 43)
    assert command.shape == power.shape == energy.shape == (9001,)
    np.testing.assert_allclose(command_time, np.arange(9001) * .025,
                               rtol=0, atol=1e-9)
    np.testing.assert_allclose(steps, [[.05, .05, .025, .025]],
                               rtol=0, atol=1e-12)
    np.testing.assert_array_equal(components, baseline["components"])
    np.testing.assert_array_equal(wave, baseline["wave"])
    np.testing.assert_allclose(body[:200],
                               short["25ms_WECCCOMP_Nonlinear_Model_Predictive_body1"],
                               rtol=0, atol=1e-9)
    np.testing.assert_allclose(pto_raw[:399],
                               short["25ms_WECCCOMP_Nonlinear_Model_Predictive_pto1"],
                               rtol=0, atol=1e-9)

    # Keep the archived 50 ms controller cadence as a negative control.
    half_step_change = np.abs(command[1::2] - command[:-1:2])
    after_30s = command_time[1::2] >= 30
    assert np.count_nonzero(half_step_change[after_30s] > 1e-6) > .99 * 3900
    assert _difference(body[:200, 5], baseline["body"][:200, 5])["maximum"] > .0017

    # BEMIO regeneration may change floating-point fit entries at 1e-12;
    # compare every numerical dataset rather than the container hash.
    with h5py.File(HYDRO) as generated, \
            h5py.File(ROOT / "fixtures/matlab/wavestar_fixed_input_step/wavestar.h5") as pinned:
        generated_keys, pinned_keys = [], []
        generated.visititems(lambda key, obj: generated_keys.append(key)
                             if isinstance(obj, h5py.Dataset) else None)
        pinned.visititems(lambda key, obj: pinned_keys.append(key)
                          if isinstance(obj, h5py.Dataset) else None)
        assert generated_keys == pinned_keys
        for key in pinned_keys:
            if pinned[key].dtype.kind in "fiu":
                np.testing.assert_allclose(generated[key][()], pinned[key][()],
                                           rtol=0, atol=5e-12)
            else:
                np.testing.assert_array_equal(generated[key][()], pinned[key][()])

    metrics = {}
    for window, count in (("first_10s", 200), ("first_30s", 601),
                          ("full_225s", 4501)):
        for name, old, new in (
            ("pitch_rad", baseline["body"][:, 5], body[:, 5]),
            ("pitch_speed_rad_s", baseline["body"][:, 11], body[:, 11]),
            ("pto_force_N", baseline["pto"][:, 33], pto[:, 33]),
            ("command_Nm", baseline["command"], command[::2]),
        ):
            metrics[f"{window}_{name}"] = _difference(old[:count], new[:count])
    old_work = np.trapezoid(baseline["pto"][:, 33] * baseline["pto"][:, 9],
                            x=baseline["pto"][:, 0])
    new_work = np.trapezoid(pto_raw[:, 33] * pto_raw[:, 9], x=pto_raw[:, 0])
    metrics["relative_pto_work_change"] = float(
        abs(new_work - old_work) / max(abs(old_work), 1e-12))
    metrics["relative_reported_energy_change"] = float(
        abs(energy[-1] - baseline["energy"][-1])
        / max(abs(baseline["energy"][-1]), 1e-12))
    print("Full WaveStar source self-change:", json.dumps(metrics, sort_keys=True))

"""Verify the derived fine-plant, published-rate WaveStar source run."""

import os
from pathlib import Path

import numpy as np
import pytest
from scipy.io import loadmat

from wecsim import WaveStarNmpcActuator, WaveStarNmpcObserver


REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
pytestmark = pytest.mark.skipif(
    not REFERENCE, reason="split-step WaveStar MATLAB trace is absent",
)
PREFIX = "WECCCOMP_NMPC_SPLIT_STEP_DIAG_"


def _trace(data, name):
    signal = data[name]
    return np.asarray(signal["time"]).reshape(-1), np.asarray(
        signal["signals"]["values"]
    )


def test_split_step_source_rates_and_precontrol_identity():
    source = Path(REFERENCE)
    body = np.loadtxt(
        source / f"{PREFIX}WECCCOMP_Nonlinear_Model_Predictive_body1.csv",
        delimiter=",", ndmin=2,
    )
    pto = np.loadtxt(
        source / f"{PREFIX}WECCCOMP_Nonlinear_Model_Predictive_pto1.csv",
        delimiter=",", ndmin=2,
    )
    components = np.loadtxt(source / f"{PREFIX}components.csv", delimiter=",")
    prior_components = np.loadtxt(
        source / "WECCCOMP_NMPC_FINE_ACTIVE_DIAG_components.csv", delimiter=",",
    )
    assert body.shape == (15201, 25)
    assert pto.shape == (15201, 43)
    assert components.shape == (500, 4)
    np.testing.assert_allclose(components, prior_components, rtol=0, atol=1e-12)
    for values in (body, pto):
        assert np.isfinite(values).all()
        np.testing.assert_allclose(values[:, 0], np.arange(15201) * .001,
                                   rtol=0, atol=1e-9)

    log = loadmat(source / f"{PREFIX}controller.mat", simplify_cells=True)
    for name, width in (("cmd_ptoM", 1), ("estimated_states", 5),
                        ("AR_excM_pred", 40)):
        time, values = _trace(log, name)
        assert len(time) == 305, f"{name} ran at the plant rate"
        assert values.size == 305 * width
        assert np.isfinite(values).all()
        np.testing.assert_allclose(time, np.arange(305) * .05,
                                   rtol=0, atol=1e-9)
    command_time, command = _trace(log, "cmd_ptoM")
    command = command.reshape(-1)
    np.testing.assert_allclose(command[command_time < 10], 0.,
                               rtol=0, atol=1e-12)
    assert np.max(np.abs(command[command_time >= 10])) > .01
    _, stroke = _trace(log, "motor_displacement")
    stroke = stroke.reshape(-1)
    assert stroke.shape == (15201,)
    observer = WaveStarNmpcObserver(dt=.05)
    estimated = np.asarray(log["estimated_states"]["signals"]["values"])
    replayed_state = np.array([
        observer.step(stroke[index * 50], command[index - 1] if index else 0.)
        for index in range(305)
    ])
    np.testing.assert_allclose(replayed_state, estimated, rtol=0, atol=1e-9)
    actuator = WaveStarNmpcActuator()
    replayed_force = np.array([
        actuator.step(command[index // 50], point)
        if index % 50 == 0 else actuator.force_at(point)
        for index, point in enumerate(pto[:, 3])
    ])
    np.testing.assert_allclose(replayed_force, pto[:, 33],
                               rtol=0, atol=1e-6)
    print("split-step source: plant samples", len(body),
          "controller samples", len(command),
          "active command peak", np.max(np.abs(command[command_time >= 15])))

"""Verify the derived fine-plant, published-rate WaveStar source run."""

import os
from pathlib import Path

import numpy as np
import pytest
from scipy.io import loadmat

from wecsim import (
    WaveStarNmpcActuator, WaveStarNmpcController, WaveStarNmpcObserver,
    WaveStarNmpcPredictor,
)
from wecsim.irregularWave import jonswap_equal_energy_components
from wecsim.wavestar import run_wavestar_published


REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
HYDRO = os.environ.get("WEC_SIM_WAVESTAR_H5")
pytestmark = pytest.mark.skipif(
    not REFERENCE, reason="split-step WaveStar MATLAB trace is absent",
)
PREFIX = "WECCCOMP_NMPC_SPLIT_STEP_DIAG_"


def _trace(data, name):
    signal = data[name]
    return np.asarray(signal["time"]).reshape(-1), np.asarray(
        signal["signals"]["values"]
    )


def test_split_step_source_rates_and_actuator_timing():
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
    # The explicit multirate chart changes the source's execution order:
    # its torque reaches the PTO one 50 ms tick after the logged request.
    actuator = WaveStarNmpcActuator()
    delayed_force = np.array([
        actuator.step(command[max(0, index // 50 - 1)], point)
        if index % 50 == 0 else actuator.force_at(point)
        for index, point in enumerate(pto[:, 3])
    ])
    np.testing.assert_allclose(delayed_force, pto[:, 33],
                               rtol=0, atol=1e-6)
    print("split-step source: plant samples", len(body),
          "controller samples", len(command),
          "active command peak", np.max(np.abs(command[command_time >= 15])))


class _DelayedSampledPTO:
    """Mirror only the derived Simulink rate-boundary delay."""

    def __init__(self):
        self.observer = WaveStarNmpcObserver(dt=.05)
        self.predictor = WaveStarNmpcPredictor(dt=.05)
        self.controller = WaveStarNmpcController(dt=.05)
        self.actuator = WaveStarNmpcActuator()
        self.sample = 0
        self.previous_command = 0.
        self.pending_command = 0.
        self.commands = []

    def step(self, stroke, *, noise=0., dropout=False):
        assert noise == 0. and not dropout
        if self.sample % 50 == 0:
            instant = self.sample * .001
            state = self.observer.step(stroke, self.previous_command)
            forecast = self.predictor.step(state[-1])
            command = self.controller.step(state, forecast, instant)
            force = self.actuator.step(self.pending_command, stroke)
            self.pending_command = command
            self.previous_command = command
            self.commands.append(command)
        else:
            force = self.actuator.force_at(stroke)
        self.sample += 1
        return force

    def force_at_stroke(self, stroke):
        return self.actuator.force_at(stroke)


@pytest.mark.skipif(not HYDRO, reason="pinned WaveStar hydrodynamics are absent")
def test_split_step_independent_short_motion():
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
    log = loadmat(source / f"{PREFIX}controller.mat", simplify_cells=True)
    _, source_command = _trace(log, "cmd_ptoM")
    sea = jonswap_equal_energy_components(
        HYDRO, significant_height=.1042, peak_period=1.836,
        directions=np.array([0.]), spreading=np.array([1.]),
        gamma=3.3, phase=components[:, 3, None],
    )
    controller = _DelayedSampledPTO()
    response = run_wavestar_published(
        HYDRO, sea, dt=.001, end_time=15.2, ramp_time=25,
        g=9.80665, pto_controller=controller,
    )
    np.testing.assert_allclose(response.time, body[:, 0], rtol=0, atol=1e-9)
    errors = {
        "pitch_rad": np.max(np.abs(response.angle - body[:, 5])),
        "pitch_speed_rad_s": np.max(np.abs(
            response.angular_speed - body[:, 11])),
        "stroke_m": np.max(np.abs(response.pto_stroke - pto[:, 3])),
        "axial_force_N": np.max(np.abs(response.pto_force - pto[:, 33])),
        "requested_torque_Nm": np.max(np.abs(
            np.asarray(controller.commands) - source_command.reshape(-1))),
    }
    limits = {
        "pitch_rad": 1e-4,
        "pitch_speed_rad_s": 1.5e-3,
        "stroke_m": 2e-5,
        "axial_force_N": .1,
        "requested_torque_Nm": .025,
    }
    for name, limit in limits.items():
        assert errors[name] < limit, f"{name}: {errors[name]} >= {limit}"
    print("split-step short closed-loop errors:", errors)

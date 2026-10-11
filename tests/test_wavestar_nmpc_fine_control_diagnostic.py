"""Diagnose WaveStar resistive control with a derived 1 ms MATLAB run."""

import os
from pathlib import Path

import numpy as np
import pytest
from scipy.io import loadmat

from wecsim import (
    WaveStarNmpcActuator, WaveStarNmpcController, WaveStarNmpcObserver,
)
from wecsim.irregularWave import jonswap_equal_energy_components
from wecsim.wavestar import run_wavestar_published


REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
HYDRO = os.environ.get("WEC_SIM_WAVESTAR_H5")
pytestmark = pytest.mark.skipif(
    not (REFERENCE and HYDRO),
    reason="derived WaveStar MATLAB output and hydrodynamics are absent",
)
PREFIX = "WECCCOMP_NMPC_FINE_CTRL_DIAG_"


def _source(name):
    return np.loadtxt(Path(REFERENCE) / f"{PREFIX}{name}.csv",
                      delimiter=",", ndmin=2)


class _ResistivePto:
    def __init__(self):
        self.observer = WaveStarNmpcObserver(dt=.001)
        self.controller = WaveStarNmpcController(dt=.001)
        self.actuator = WaveStarNmpcActuator()
        self.previous_command = 0.
        self.sample = 0
        self.commands = []

    def step(self, stroke, noise=0., dropout=False):
        instant = self.sample * .001
        state = self.observer.step(stroke, self.previous_command)
        command = self.controller.step(state, np.zeros(40), instant)
        force = self.actuator.step(command, stroke)
        self.previous_command = command
        self.commands.append(command)
        self.sample += 1
        return force

    def force_at_stroke(self, stroke):
        return self.actuator.force_at(stroke)


def test_fine_step_resistive_control_diagnostic():
    components = _source("components")
    source_float = _source("WECCCOMP_Nonlinear_Model_Predictive_body1")
    source_pto = _source("WECCCOMP_Nonlinear_Model_Predictive_pto1")
    controller = loadmat(Path(REFERENCE) / f"{PREFIX}controller.mat",
                         simplify_cells=True)
    time = np.asarray(controller["cmd_ptoM"]["time"])
    command = np.asarray(controller["cmd_ptoM"]["signals"]["values"])
    stroke = np.asarray(controller["motor_displacement"]["signals"]["values"])
    estimated = np.asarray(controller["estimated_states"]["signals"]["values"])
    assert components.shape == (500, 4)
    assert time.shape == command.shape == stroke.shape == (14951,)
    assert estimated.shape == (14951, 5)
    assert source_float.shape == (14951, 25)
    assert source_pto.shape == (14951, 43)
    np.testing.assert_allclose(time, np.arange(14951) * .001,
                               rtol=0, atol=1e-9)
    np.testing.assert_allclose(source_float[:, 0], time,
                               rtol=0, atol=1e-9)
    np.testing.assert_allclose(command[time < 10], 0., rtol=0, atol=1e-12)
    np.testing.assert_allclose(command[time >= 10],
                               np.clip(-19.4 * estimated[time >= 10, 1],
                                       -12, 12),
                               rtol=0, atol=1e-11)
    observer = WaveStarNmpcObserver(dt=.001)
    replay_state = np.array([
        observer.step(position, command[index - 1] if index else 0.)
        for index, position in enumerate(stroke)
    ])
    np.testing.assert_allclose(replay_state, estimated,
                               rtol=0, atol=1e-10)
    actuator = WaveStarNmpcActuator()
    replay_force = np.array([
        actuator.step(request, position)
        for request, position in zip(command, stroke)
    ])
    np.testing.assert_allclose(replay_force, source_pto[:, 33],
                               rtol=0, atol=1e-8)

    sea = jonswap_equal_energy_components(
        HYDRO, significant_height=.1042, peak_period=1.836,
        directions=np.array([0.]), spreading=np.array([1.]),
        gamma=3.3, phase=components[:, 3, None],
    )
    pto = _ResistivePto()
    response = run_wavestar_published(
        HYDRO, sea, dt=.001, end_time=14.95, ramp_time=25,
        g=9.80665, pto_controller=pto,
    )
    assert response.time.shape == time.shape
    precontrol = time < 10
    assert np.max(np.abs(
        response.angle[precontrol] - source_float[precontrol, 5]
    )) < 2e-5
    source_forces = _source("body1_forces")
    assert source_forces.shape == (14951, 37)
    errors = {
        "pitch_rad": np.max(np.abs(response.angle - source_float[:, 5])),
        "pitch_speed_rad_s": np.max(np.abs(
            response.angular_speed - source_float[:, 11])),
        "float_center_m": np.max(np.abs(
            response.float_position[:, [0, 2]] - source_float[:, [1, 3]])),
        "float_center_speed_m_s": np.max(np.abs(
            response.float_velocity[:, [0, 2]] - source_float[:, [7, 9]])),
        "pto_stroke_m": np.max(np.abs(
            response.pto_stroke - source_pto[:, 3])),
        "pto_speed_m_s": np.max(np.abs(
            response.pto_speed - source_pto[:, 9])),
        "command_Nm": np.max(np.abs(np.asarray(pto.commands) - command)),
        "pto_force_N": np.max(np.abs(
            response.pto_force - source_pto[:, 33])),
        "radiation_N_or_Nm": np.max(np.abs(
            response.radiation_force - source_forces[:, 1:7])),
    }
    limits = {
        "pitch_rad": 5e-5,
        "pitch_speed_rad_s": 4e-4,
        "float_center_m": 2.5e-5,
        "float_center_speed_m_s": 2.5e-4,
        "pto_stroke_m": 1e-5,
        "pto_speed_m_s": 1e-4,
        "command_Nm": 5e-3,
        "pto_force_N": 3e-2,
        "radiation_N_or_Nm": 1e-3,
    }
    for name, limit in limits.items():
        assert errors[name] < limit, f"{name}: {errors[name]} >= {limit}"
    print("fine-step resistive-control errors:", errors)

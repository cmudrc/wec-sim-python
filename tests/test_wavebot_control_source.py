"""Pinned MATLAB WaveBot controller inputs and sampled-force laws."""

import os
from pathlib import Path

import numpy as np
import pytest
from scipy.io import loadmat

from wecsim.environment.irregularWave import (
    pm_equal_energy_components, synthesize_irregular_response,
)
from wecsim.control.wavebotControl import (
    WaveBotPIController, WaveBotPowerAutotuner, run_wavebot_power_control,
)


REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
HYDRO = os.environ.get("WEC_SIM_WAVEBOT_H5")
APPLICATION = os.environ.get("WEC_SIM_WAVEBOT_CONTROL_DIR")
pytestmark = pytest.mark.skipif(
    not (REFERENCE and HYDRO),
    reason="pinned WaveBot controller output and hydrodynamics are absent",
)


def _source(name):
    values = np.loadtxt(
        Path(REFERENCE) / f"WAVEBOT_CONTROL_{name}.csv",
        delimiter=",", ndmin=2,
    )
    assert np.isfinite(values).all(), name
    return values


def _components():
    phase = _source("components")[:, 3, None]
    return pm_equal_energy_components(
        HYDRO, significant_height=.254, peak_period=3.5,
        directions=np.array([0.]), spreading=np.array([1.]),
        phase=phase,
    )


def _tuner_inputs():
    if not APPLICATION:
        pytest.skip("pinned WaveBot controller impedance is absent")
    base = Path(APPLICATION)
    return (loadmat(base / "Z_included.mat")["Z"],
            loadmat(base / "f_vec.mat")["f_vec"].ravel())


def test_published_pm_sea_and_excitation():
    source_components = _source("components")
    source_wave = _source("wave")
    source_force = _source("body1_forces")
    assert source_components.shape == (500, 4)
    assert source_wave.shape == (10_001, 2)
    assert source_force.shape == (10_001, 43)

    components = _components()
    np.testing.assert_allclose(components.omega, source_components[:, 0],
                               rtol=0, atol=1e-12)
    np.testing.assert_allclose(components.spectral_amplitude,
                               source_components[:, 1], rtol=0, atol=1e-12)
    np.testing.assert_allclose(components.d_omega, source_components[:, 2],
                               rtol=0, atol=1e-12)
    incident = synthesize_irregular_response(
        HYDRO, components, dt=.01, end_time=100, ramp_time=20,
        rho=1025, g=9.81,
    )
    np.testing.assert_allclose(incident.time, source_wave[:, 0],
                               rtol=0, atol=1e-10)
    np.testing.assert_allclose(incident.elevation, source_wave[:, 1],
                               rtol=0, atol=1e-10)
    np.testing.assert_allclose(incident.excitation_force,
                               source_force[:, 1:7], rtol=0, atol=1e-7)


def test_published_mooring_drag_and_indexed_damping():
    body = _source("ControlTests_body1")
    force = _source("body1_forces")
    mooring = _source("mooring")
    assert body.shape == (10_001, 25)
    assert mooring.shape == (10_001, 19)
    np.testing.assert_allclose(body[:, 0], np.arange(10_001) * .01,
                               rtol=0, atol=1e-8)
    np.testing.assert_allclose(mooring[:, 0], body[:, 0],
                               rtol=0, atol=1e-8)
    stiffness = np.diag([24_000, 24_000, 5_000, 0, 1_000, 0])
    damping = np.diag([5, 0, 5, 0, 0, 0])
    expected_mooring = (-mooring[:, 1:7] @ stiffness.T
                        - mooring[:, 7:13] @ damping.T)
    np.testing.assert_allclose(mooring[:, 13:19], expected_mooring,
                               rtol=0, atol=1e-6)

    velocity = body[:, 7:13]
    drag = (.5 * 1025 * np.array([1.15, 1.15, 1, .5, .5, 0])
            * np.array([2.9568, 2.9568, 5.4739, 5.4739, 5.4739, 0]))
    np.testing.assert_allclose(force[:, 25:31],
                               drag * np.abs(velocity) * velocity,
                               rtol=0, atol=1e-6)
    indexed_damping = np.zeros((6, 6))
    indexed_damping[[0, 2, 4], 0] = [1000, 1000, 100]
    np.testing.assert_allclose(force[:, 31:37],
                               velocity @ indexed_damping.T,
                               rtol=0, atol=1e-6)


def test_sampled_pi_force_reconstructs_logged_actuation():
    body = _source("ControlTests_body1")
    gains = _source("gains_real")
    actuation = _source("F_actuation_real")
    excitation = (_source("FexcOut_real")[:, 1:]
                  + 1j * _source("FexcOut_imag")[:, 1:])
    weights = _source("WEIGHTS_real")
    assert gains.shape == (10_001, 7)
    assert actuation.shape[1] == 4
    assert excitation.shape == (13, 615)
    assert weights.shape == (13, 10)
    np.testing.assert_allclose(weights[:, 1:4], 25, rtol=0, atol=1e-12)
    np.testing.assert_allclose(weights[:, 4:], 0, rtol=0, atol=1e-12)
    np.testing.assert_allclose(weights[:, 0], np.arange(13) * 8,
                               rtol=0, atol=1e-8)
    assert np.max(np.abs(excitation[1:])) > 10

    velocity = body[:, [7, 9, 11]]
    controller = WaveBotPIController()
    predicted = np.empty((len(body), 3))
    for index in range(len(body)):
        controller.set_gains(gains[index, 1:])
        predicted[index] = controller.force(velocity[index])
        controller.advance(velocity[index], .01)
    logged = np.column_stack([
        np.interp(body[:, 0], actuation[:, 0], actuation[:, channel])
        for channel in (1, 2, 3)
    ])
    np.testing.assert_allclose(predicted, logged, rtol=0, atol=1e-8)


def test_published_spectral_estimator_at_all_update_times():
    impedance, frequency = _tuner_inputs()
    tuner = WaveBotPowerAutotuner(impedance, frequency)
    body = _source("ControlTests_body1")
    actuation = _source("F_actuation_real")
    target = (_source("FexcOut_real")[:, 1:]
              + 1j * _source("FexcOut_imag")[:, 1:])
    assert target.shape == (13, 615)
    for update in range(13):
        sample_time = np.arange(0, update * 8, .25)
        velocity = np.column_stack([
            np.interp(sample_time, body[:, 0], body[:, column])
            for column in (7, 9, 11)
        ])
        delayed = np.column_stack([
            np.interp(np.maximum(sample_time - .01, 0),
                      actuation[:, 0], actuation[:, column])
            for column in (1, 2, 3)
        ])
        observed = tuner.excitation_spectrum(velocity, delayed)
        expected = target[update].reshape(205, 3, order="F")
        np.testing.assert_allclose(observed, expected, rtol=0, atol=2e-12)


def test_independent_adaptive_power_control_trajectory():
    impedance, frequency = _tuner_inputs()
    response = run_wavebot_power_control(
        HYDRO, _components(), impedance, frequency,
        source_linear_damping=True,
    )
    body = _source("ControlTests_body1")
    gains = _source("gains_real")
    actuation = _source("F_actuation_real")
    excitation = (_source("FexcOut_real")[:, 1:]
                  + 1j * _source("FexcOut_imag")[:, 1:])
    np.testing.assert_allclose(response.time, body[:, 0], rtol=0, atol=1e-8)
    np.testing.assert_allclose(response.gains[0], 0, rtol=0, atol=1e-12)
    np.testing.assert_allclose(response.gains[1],
                               [-1000, 0, -1000, 500, -200, 0],
                               rtol=0, atol=1e-12)
    for coordinate, position, velocity in ((0, 1, 7), (2, 3, 9), (4, 5, 11)):
        np.testing.assert_allclose(response.body_position[:, coordinate],
                                   body[:, position], rtol=0, atol=5e-4)
        np.testing.assert_allclose(response.body_velocity[:, coordinate],
                                   body[:, velocity], rtol=0, atol=2e-3)
    logged_command = np.column_stack([
        np.interp(response.time, actuation[:, 0], actuation[:, channel])
        for channel in (1, 2, 3)
    ])
    np.testing.assert_allclose(response.actuator_command, logged_command,
                               rtol=0, atol=2.)
    np.testing.assert_allclose(response.gains, gains[:, 1:],
                               rtol=0, atol=40.)
    assert response.excitation_estimates.shape == (13, 205, 3)
    expected = np.stack([row.reshape(205, 3, order="F")
                         for row in excitation])
    np.testing.assert_allclose(response.excitation_estimates, expected,
                               rtol=0, atol=.07)

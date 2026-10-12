"""Published WaveStar fault application's source forces and linkage."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim.environment.irregularWave import (
    jonswap_equal_energy_components, synthesize_irregular_response,
)
from wecsim.dynamics.friction import StribeckFriction
from wecsim.dynamics.wavestar import WaveStarLinkage
from wecsim.control.wavestarFault import WaveStarFaultController


REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
HYDRO = os.environ.get("WEC_SIM_WAVESTAR_H5")
pytestmark = pytest.mark.skipif(
    not (REFERENCE and HYDRO),
    reason="pinned WaveStar fault output and hydrodynamics are absent",
)


def _source(name):
    values = np.loadtxt(
        Path(REFERENCE) / f"WECCCOMP_FAULT_SOURCE_{name}.csv",
        delimiter=",", ndmin=2,
    )
    assert np.isfinite(values).all()
    return values


def test_fault_application_wave_and_excitation():
    source_components = _source("components")
    source_wave = _source("wave")
    source_float = _source("WECCCOMP_Fault_Implementation_body1")
    assert source_components.shape == (500, 4)
    assert source_wave.shape == (14_121, 2)
    components = jonswap_equal_energy_components(
        HYDRO, significant_height=.0625, peak_period=1.412,
        directions=np.array([0.]), spreading=np.array([1.]),
        phase=source_components[:, 3, None],
    )
    np.testing.assert_allclose(components.omega, source_components[:, 0],
                               rtol=0, atol=1e-10)
    np.testing.assert_allclose(components.spectral_amplitude,
                               source_components[:, 1], rtol=0, atol=1e-12)
    np.testing.assert_allclose(components.d_omega, source_components[:, 2],
                               rtol=0, atol=1e-10)
    incident = synthesize_irregular_response(
        HYDRO, components, dt=.01, end_time=141.2,
        ramp_time=7.06, rho=1000, g=9.81,
    )
    np.testing.assert_allclose(incident.elevation, source_wave[:, 1],
                               rtol=0, atol=1e-10)
    np.testing.assert_allclose(incident.excitation_force,
                               source_float[:, 19:25], rtol=0, atol=1e-7)


def test_fault_application_linkage_and_friction_window():
    float_body = _source("WECCCOMP_Fault_Implementation_body1")
    arm = _source("WECCCOMP_Fault_Implementation_body2")
    frame = _source("WECCCOMP_Fault_Implementation_body3")
    pto = [_source(f"WECCCOMP_Fault_Implementation_pto{index}")
           for index in range(1, 5)]
    assert all(values.shape == (14_121, 25)
               for values in (float_body, arm, frame))
    assert all(values.shape == (14_121, 43) for values in pto)
    np.testing.assert_allclose(float_body[:, 0], np.arange(14_121) * .01,
                               rtol=0, atol=1e-8)
    np.testing.assert_allclose(arm[:, 5], float_body[:, 5],
                               rtol=0, atol=1e-12)
    assert np.max(np.abs(frame[:, 1:13] - frame[0, 1:13])) < 1e-12
    linkage = WaveStarLinkage()
    angle = float_body[:, 5]
    speed = float_body[:, 11]
    for response in (float_body, arm):
        neutral = response[0, [1, 3]]
        np.testing.assert_allclose(
            linkage.point_position(neutral, angle),
            response[:, [1, 3]], rtol=0, atol=1e-12,
        )
        np.testing.assert_allclose(
            linkage.point_velocity(neutral, angle, speed),
            response[:, [7, 9]], rtol=0, atol=1e-12,
        )
    np.testing.assert_allclose(linkage.pto_stroke(angle), pto[0][:, 3],
                               rtol=0, atol=1e-12)
    np.testing.assert_allclose(linkage.pto_speed(angle, speed), pto[0][:, 9],
                               rtol=0, atol=1e-12)
    np.testing.assert_allclose(pto[3][:, 5], angle, rtol=0, atol=1e-12)
    friction = StribeckFriction(
        breakaway_torque=.25, breakaway_speed=.1,
        coulomb_torque=.2, viscous_damping=.001,
    )
    for joint in pto[1:3]:
        np.testing.assert_allclose(
            friction.torque(joint[:, 11]), joint[:, 35],
            rtol=0, atol=1e-12,
        )
    time = pto[3][:, 0]
    pivot_torque = pto[3][:, 35]
    active = (time > 55) & (time < 115)
    assert np.max(np.abs(pivot_torque[active])) > .2
    np.testing.assert_allclose(
        friction.torque(pto[3][active, 11]), pivot_torque[active],
        rtol=0, atol=1e-12,
    )
    np.testing.assert_allclose(
        pivot_torque[~active],
        0, rtol=0, atol=1e-12,
    )
    assert np.max(np.abs(pto[0][:, 33])) > 1


def test_fault_sensor_and_sampled_controller():
    true = _source("fault_true_position")
    measured = _source("fault_measured_position")
    motor = _source("WECCCOMP_Fault_Implementation_pto1")
    assert true.shape == measured.shape == (141_201, 2)
    np.testing.assert_allclose(true[:, 0], np.arange(141_201) * .001,
                               rtol=0, atol=1e-8)
    np.testing.assert_allclose(measured[:, 0], true[:, 0],
                               rtol=0, atol=1e-12)
    dropout = measured[:, 1] == 0
    assert .02 < np.mean(dropout) < .04
    noise = measured[:, 1] - true[:, 1]
    assert .0025 < np.std(noise[~dropout]) < .0035

    angle = true[:, 1] + 1.170165
    length = np.sqrt(.412**2 + .2**2 - 2 * .412 * .2 * np.cos(angle))
    stroke = length - .381408
    np.testing.assert_allclose(
        WaveStarFaultController.true_angle(motor[:, 3]), true[::10, 1],
        rtol=0, atol=1e-12,
    )
    controller = WaveStarFaultController()
    force = np.array([
        controller.step(stroke[index], noise=noise[index],
                        dropout=dropout[index])
        for index in range(len(stroke))
    ])
    np.testing.assert_allclose(force[::10], motor[:, 33],
                               rtol=0, atol=1e-7)

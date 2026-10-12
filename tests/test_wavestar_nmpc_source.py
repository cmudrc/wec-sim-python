"""Check the pinned WaveStar NMPC application's sea and mechanical signals."""

import os
from pathlib import Path

import h5py
import numpy as np
import pytest
from scipy.io import loadmat

from wecsim.environment.irregularWave import (
    jonswap_equal_energy_components, synthesize_irregular_response,
)
from wecsim.dynamics.wavestar import WaveStarLinkage, run_wavestar_published
from wecsim import (
    WaveStarNmpcActuator, WaveStarNmpcController,
    WaveStarNmpcObserver, WaveStarNmpcPredictor, WaveStarNmpcPTO,
)


REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
HYDRO = os.environ.get("WEC_SIM_WAVESTAR_H5")
pytestmark = pytest.mark.skipif(
    not (REFERENCE and HYDRO),
    reason="pinned WaveStar NMPC output and hydrodynamics are absent",
)


def _source(name):
    values = np.loadtxt(
        Path(REFERENCE) / f"WECCCOMP_NMPC_SOURCE_{name}.csv",
        delimiter=",", ndmin=2,
    )
    assert np.isfinite(values).all()
    return values


def test_nmpc_source_radiation_fit_at_observed_high_pitch():
    """The source fit can add low-frequency energy at its own large pitch."""
    source_float = _source("WECCCOMP_Nonlinear_Model_Predictive_body1")
    assert source_float.shape == (4501, 25)
    high_pitch = float(np.max(source_float[:, 5]))
    assert .35 < high_pitch < .4
    high_samples = source_float[:, 5] > .25
    assert 0 < np.count_nonzero(high_samples) < .01 * len(source_float)
    assert source_float[np.flatnonzero(high_samples)[0], 0] > 15

    with h5py.File(HYDRO) as h5:
        omega = float(np.asarray(h5["simulation_parameters/w"]).ravel()[0])
        rho = float(np.asarray(h5["simulation_parameters/rho"]).item())
        center = np.asarray(h5["body1/properties/cg"]).ravel()
        radiation = h5["body1/hydro_coeffs/radiation_damping"]
        bem = np.asarray(radiation["all"])[:, :, 0] * rho
        fit = radiation["state_space"]
        order = np.asarray(fit["it"], dtype=int)
        matrix_a = np.asarray(fit["A/all"])
        matrix_b = np.asarray(fit["B/all"])
        matrix_c = np.asarray(fit["C/all"]) * rho
    assert .15 < omega < .25

    def projected_damping(angle):
        jacobian = np.zeros(6)
        jacobian[[0, 2]] = WaveStarLinkage().point_velocity(
            center[[0, 2]], angle, 1.,
        )
        jacobian[4] = 1
        fitted = 0j
        for output in (0, 2, 4):
            for input_axis in (0, 2, 4):
                count = order[output, input_axis]
                if count:
                    state_a = matrix_a[output, input_axis, :count, :count]
                    state_b = matrix_b[output, input_axis, :count, 0]
                    state_c = matrix_c[output, input_axis, 0, :count]
                    fitted += (jacobian[output] * jacobian[input_axis]
                               * state_c @ np.linalg.solve(
                                   1j * omega * np.eye(count) - state_a,
                                   state_b,
                               ))
        return fitted.real, float(jacobian @ bem @ jacobian)

    neutral_fit, _ = projected_damping(0.)
    high_fit, high_bem = projected_damping(high_pitch)
    assert neutral_fit > .1
    assert high_fit < -.05
    assert high_bem > .05


def test_nmpc_published_sea_and_excitation():
    saved = _source("components")
    wave = _source("wave")
    source_float = _source("WECCCOMP_Nonlinear_Model_Predictive_body1")
    assert saved.shape == (500, 4)
    assert wave.shape == (4501, 2)
    assert source_float.shape[1] == 25
    np.testing.assert_allclose(wave[:, 0], np.arange(4501) * .05,
                               rtol=0, atol=1e-9)
    components = jonswap_equal_energy_components(
        HYDRO, significant_height=.1042, peak_period=1.836,
        directions=np.array([0.]), spreading=np.array([1.]),
        gamma=3.3, phase=saved[:, 3, None],
    )
    np.testing.assert_allclose(components.omega, saved[:, 0],
                               rtol=0, atol=1e-10)
    np.testing.assert_allclose(components.spectral_amplitude, saved[:, 1],
                               rtol=0, atol=1e-11)
    np.testing.assert_allclose(components.d_omega, saved[:, 2],
                               rtol=0, atol=1e-10)
    incident = synthesize_irregular_response(
        HYDRO, components, dt=.05, end_time=225, ramp_time=25,
        rho=1000, g=9.80665,
    )
    np.testing.assert_allclose(incident.elevation, wave[:, 1],
                               rtol=0, atol=1e-10)
    indices = np.rint(source_float[:, 0] / .05).astype(int)
    np.testing.assert_allclose(source_float[:, 0], wave[indices, 0],
                               rtol=0, atol=1e-9)
    np.testing.assert_allclose(
        incident.excitation_force[indices], source_float[:, 19:25],
        rtol=0, atol=1e-7,
    )


def test_nmpc_published_linkage_and_control_trace():
    source_float = _source("WECCCOMP_Nonlinear_Model_Predictive_body1")
    arm = _source("WECCCOMP_Nonlinear_Model_Predictive_body2")
    pto = _source("WECCCOMP_Nonlinear_Model_Predictive_pto1")
    assert arm.shape == source_float.shape
    assert pto.shape[1] == 43
    np.testing.assert_allclose(arm[:, 0], source_float[:, 0],
                               rtol=0, atol=1e-9)
    linkage = WaveStarLinkage()
    angle = source_float[:, 5]
    speed = source_float[:, 11]
    for source in (source_float, arm):
        neutral = source[0, [1, 3]]
        np.testing.assert_allclose(
            linkage.point_position(neutral, angle),
            source[:, [1, 3]], rtol=0, atol=1e-10,
        )
        np.testing.assert_allclose(
            linkage.point_velocity(neutral, angle, speed),
            source[:, [7, 9]], rtol=0, atol=1e-9,
        )
    np.testing.assert_allclose(linkage.pto_stroke(angle), pto[:, 3],
                               rtol=0, atol=1e-9)
    np.testing.assert_allclose(linkage.pto_speed(angle, speed), pto[:, 9],
                               rtol=0, atol=1e-9)
    assert np.max(np.abs(pto[pto[:, 0] > 15, 33])) > 1

    controller = loadmat(
        Path(REFERENCE) / "WECCCOMP_NMPC_SOURCE_controller.mat",
        simplify_cells=True,
    )
    for name in ("cmd_ptoM", "Output_power", "Output_energy",
                 "estimated_states", "AR_excM_pred", "excM_wecSim",
                 "motor_displacement"):
        assert name in controller, f"source controller trace {name} missing"


def test_nmpc_actuator_on_source_command_and_stroke():
    pto = _source("WECCCOMP_Nonlinear_Model_Predictive_pto1")
    controller = loadmat(
        Path(REFERENCE) / "WECCCOMP_NMPC_SOURCE_controller.mat",
        simplify_cells=True,
    )
    command = np.asarray(controller["cmd_ptoM"]["signals"]["values"])
    stroke = np.asarray(
        controller["motor_displacement"]["signals"]["values"]
    )
    time = np.asarray(controller["cmd_ptoM"]["time"])
    assert command.shape == stroke.shape == time.shape == (4501,)
    np.testing.assert_allclose(time, np.arange(4501) * .05,
                               rtol=0, atol=1e-9)
    assert np.max(np.abs(command[time < 10])) < 1e-12
    assert np.max(np.abs(command)) <= 12 + 1e-10
    assert np.max(np.abs(command[time >= 15])) > 1
    estimated = np.asarray(
        controller["estimated_states"]["signals"]["values"]
    )
    assert estimated.shape == (4501, 5)
    resistive = (time >= 10) & (time < 15)
    np.testing.assert_allclose(
        command[resistive], np.clip(-19.4 * estimated[resistive, 1], -12, 12),
        rtol=0, atol=1e-12,
    )
    actuator = WaveStarNmpcActuator()
    force = np.array([
        actuator.step(request, position)
        for request, position in zip(command, stroke)
    ])
    indices = np.rint(pto[:, 0] / .05).astype(int)
    np.testing.assert_allclose(stroke[indices], pto[:, 3],
                               rtol=0, atol=1e-12)
    np.testing.assert_allclose(force[indices], pto[:, 33],
                               rtol=0, atol=1e-7)


def test_nmpc_observer_on_source_command_and_stroke():
    controller = loadmat(
        Path(REFERENCE) / "WECCCOMP_NMPC_SOURCE_controller.mat",
        simplify_cells=True,
    )
    command = np.asarray(controller["cmd_ptoM"]["signals"]["values"])
    stroke = np.asarray(
        controller["motor_displacement"]["signals"]["values"]
    )
    saved = np.asarray(
        controller["estimated_states"]["signals"]["values"]
    )
    assert command.shape == stroke.shape == (4501,)
    assert saved.shape == (4501, 5)
    observer = WaveStarNmpcObserver()
    replay = np.array([
        observer.step(position, command[index - 1] if index else 0.)
        for index, position in enumerate(stroke)
    ])
    np.testing.assert_allclose(replay, saved, rtol=0, atol=1e-10)


def test_nmpc_ar_forecast_from_estimated_excitation():
    controller = loadmat(
        Path(REFERENCE) / "WECCCOMP_NMPC_SOURCE_controller.mat",
        simplify_cells=True,
    )
    moment = np.asarray(
        controller["estimated_states"]["signals"]["values"]
    )[:, 4]
    source_forecast = np.asarray(
        controller["AR_excM_pred"]["signals"]["values"]
    )
    assert moment.shape == (4501,)
    assert source_forecast.shape == (40, 4501)
    predictor = WaveStarNmpcPredictor()
    replay = np.column_stack([predictor.step(value) for value in moment])
    np.testing.assert_allclose(replay[:, :200], source_forecast[:, :200],
                               rtol=0, atol=1e-12)
    # The first fit produces an unstable 10.05-10.95 s source forecast
    # exceeding 5e6 N m. Its long horizon amplifies tiny solver differences.
    np.testing.assert_allclose(replay[:, 200:220],
                               source_forecast[:, 200:220],
                               rtol=0, atol=.25)
    np.testing.assert_allclose(replay[:, 220:],
                               source_forecast[:, 220:],
                               rtol=0, atol=1e-9)

    # Compose both Python components from measured source inputs. After the
    # first AR refit at 11 s, the forecast remains tightly paired without
    # feeding either saved estimator state or saved forecast into Python.
    command = np.asarray(controller["cmd_ptoM"]["signals"]["values"])
    stroke = np.asarray(
        controller["motor_displacement"]["signals"]["values"]
    )
    observer = WaveStarNmpcObserver()
    predictor = WaveStarNmpcPredictor()
    composed = np.column_stack([
        predictor.step(observer.step(
            position, command[index - 1] if index else 0.,
        )[-1])
        for index, position in enumerate(stroke)
    ])
    np.testing.assert_allclose(composed[:, 220:],
                               source_forecast[:, 220:],
                               rtol=0, atol=1e-9)
    assert np.max(np.abs(composed[:, 200:220]
                         - source_forecast[:, 200:220])) < 25


def test_nmpc_controller_command_from_source_and_composed_inputs():
    controller = loadmat(
        Path(REFERENCE) / "WECCCOMP_NMPC_SOURCE_controller.mat",
        simplify_cells=True,
    )
    time = np.asarray(controller["cmd_ptoM"]["time"])
    source_command = np.asarray(
        controller["cmd_ptoM"]["signals"]["values"]
    )
    source_state = np.asarray(
        controller["estimated_states"]["signals"]["values"]
    )
    source_forecast = np.asarray(
        controller["AR_excM_pred"]["signals"]["values"]
    )
    stroke = np.asarray(
        controller["motor_displacement"]["signals"]["values"]
    )
    assert time.shape == source_command.shape == stroke.shape == (4501,)
    assert source_state.shape == (4501, 5)
    assert source_forecast.shape == (40, 4501)

    control = WaveStarNmpcController()
    replay = np.array([
        control.step(source_state[index], source_forecast[:, index], instant)
        for index, instant in enumerate(time)
    ])
    # Two source exports with identical wave components take slightly
    # different RTI nominal-horizon paths. Both retain a tight aggregate
    # command match; the larger pointwise error is near a zero crossing.
    replay_error = np.abs(replay - source_command)
    assert np.max(replay_error) < .01
    assert np.sqrt(np.mean(replay_error**2)) < .0003

    observer = WaveStarNmpcObserver()
    predictor = WaveStarNmpcPredictor()
    control = WaveStarNmpcController()
    composed = []
    for index, (instant, position) in enumerate(zip(time, stroke)):
        state = observer.step(position, source_command[index - 1]
                              if index else 0.)
        forecast = predictor.step(state[-1])
        composed.append(control.step(state, forecast, instant))
    composed_error = np.abs(np.asarray(composed) - source_command)
    assert np.max(composed_error) < .01
    assert np.sqrt(np.mean(composed_error**2)) < .0003


def test_nmpc_unforced_plant_motion_diagnostic():
    """Bound the pre-control solver difference without claiming NMPC parity."""
    source_float = _source("WECCCOMP_Nonlinear_Model_Predictive_body1")
    source_pto = _source("WECCCOMP_Nonlinear_Model_Predictive_pto1")
    source_forces = _source("body1_forces")
    assert source_float.shape == (4501, 25)
    assert source_pto.shape == (4501, 43)
    assert source_forces.shape == (4501, 37)
    components = _source("components")
    sea = jonswap_equal_energy_components(
        HYDRO, significant_height=.1042, peak_period=1.836,
        directions=np.array([0.]), spreading=np.array([1.]),
        gamma=3.3, phase=components[:, 3, None],
    )
    # The published source uses ode8 at 0.05 s. A finer Python step keeps
    # the fitted radiation system stable; the numerical paths still differ.
    response = run_wavestar_published(
        HYDRO, sea, dt=.001, end_time=9.95, ramp_time=25,
        g=9.80665, output_stride=50,
    )
    count = len(response.time)
    np.testing.assert_allclose(response.time, source_float[:count, 0],
                               rtol=0, atol=1e-9)
    assert np.max(np.abs(source_pto[:count, 31:37])) < 1e-12
    assert np.max(np.abs(response.angle - source_float[:count, 5])) < .004
    assert np.max(np.abs(
        response.angular_speed - source_float[:count, 11]
    )) < .04
    assert np.max(np.abs(
        response.float_position[:, [0, 2]]
        - source_float[:count][:, [1, 3]]
    )) < .002
    assert np.max(np.abs(
        response.pto_stroke - source_pto[:count, 3]
    )) < .001
    assert np.max(np.abs(
        response.radiation_force - source_forces[:count, 1:7]
    )) < .2


def test_nmpc_moving_arm_short_window_step_refinement():
    """The moving PTO arm must not create plant-step drift before 30 s."""
    components = _source("components")
    sea = jonswap_equal_energy_components(
        HYDRO, significant_height=.1042, peak_period=1.836,
        directions=np.array([0.]), spreading=np.array([1.]),
        gamma=3.3, phase=components[:, 3, None],
    )
    runs = []
    for plant_dt in (.001, .0005):
        pto = WaveStarNmpcPTO(plant_dt=plant_dt, control_dt=.05)
        response = run_wavestar_published(
            HYDRO, sea, dt=plant_dt, end_time=30, ramp_time=25,
            g=9.80665, pto_controller=pto,
            output_stride=round(.05 / plant_dt),
        )
        runs.append((response, pto))
    (coarse, coarse_pto), (fine, fine_pto) = runs
    np.testing.assert_allclose(coarse.time, fine.time, rtol=0, atol=1e-10)
    assert np.max(np.abs(
        coarse_pto.command_torque[coarse.time >= 15]
    )) > 1
    assert np.max(np.abs(coarse.angle - fine.angle)) < 1e-5
    assert np.max(np.abs(
        coarse_pto.command_torque - fine_pto.command_torque
    )) < 1e-3
    assert np.max(np.abs(coarse.pto_force - fine.pto_force)) < .01


def test_nmpc_independent_full_closed_loop_diagnostic():
    """Expose the published coarse-solver gap without tuning plant physics."""
    source_float = _source("WECCCOMP_Nonlinear_Model_Predictive_body1")
    source_pto = _source("WECCCOMP_Nonlinear_Model_Predictive_pto1")
    components = _source("components")
    controller = loadmat(
        Path(REFERENCE) / "WECCCOMP_NMPC_SOURCE_controller.mat",
        simplify_cells=True,
    )
    source_command = np.asarray(
        controller["cmd_ptoM"]["signals"]["values"]
    )
    sea = jonswap_equal_energy_components(
        HYDRO, significant_height=.1042, peak_period=1.836,
        directions=np.array([0.]), spreading=np.array([1.]),
        gamma=3.3, phase=components[:, 3, None],
    )
    pto = WaveStarNmpcPTO(plant_dt=.001, control_dt=.05)
    response = run_wavestar_published(
        HYDRO, sea, dt=.001, end_time=225, ramp_time=25,
        g=9.80665, pto_controller=pto, output_stride=50,
    )
    assert source_float.shape == (4501, 25)
    assert source_pto.shape == (4501, 43)
    assert source_command.shape == (4501,)
    np.testing.assert_allclose(response.time, source_float[:, 0],
                               rtol=0, atol=1e-9)
    np.testing.assert_allclose(pto.command_time, response.time,
                               rtol=0, atol=1e-9)
    assert pto.command_torque.shape == (4501,)
    assert np.isfinite(response.angle).all()
    assert np.isfinite(response.pto_force).all()
    assert np.isfinite(pto.command_torque).all()
    assert np.max(np.abs(pto.command_torque)) <= 12 + 1e-10
    errors = {
        "pitch_rad": np.max(np.abs(response.angle - source_float[:, 5])),
        "pitch_speed_rad_s": np.max(np.abs(
            response.angular_speed - source_float[:, 11])),
        "pto_stroke_m": np.max(np.abs(
            response.pto_stroke - source_pto[:, 3])),
        "command_Nm": np.max(np.abs(
            pto.command_torque - source_command)),
        "axial_force_N": np.max(np.abs(
            response.pto_force - source_pto[:, 33])),
    }
    print("published coarse-step NMPC motion differences:", errors)


def test_nmpc_prescribed_torque_full_plant_diagnostic():
    """Isolate the plant with source torque and a moving rod-arm force."""
    source_float = _source("WECCCOMP_Nonlinear_Model_Predictive_body1")
    source_pto = _source("WECCCOMP_Nonlinear_Model_Predictive_pto1")
    components = _source("components")
    controller = loadmat(
        Path(REFERENCE) / "WECCCOMP_NMPC_SOURCE_controller.mat",
        simplify_cells=True,
    )
    command = np.asarray(controller["cmd_ptoM"]["signals"]["values"])
    assert source_float.shape == (4501, 25)
    assert source_pto.shape == (4501, 43)
    assert command.shape == (4501,)

    class PrescribedTorquePTO:
        plant_dt = .001

        def __init__(self):
            self.sample = 0
            self.actuator = WaveStarNmpcActuator()

        def step(self, stroke, *, noise=0., dropout=False):
            assert noise == 0. and not dropout
            if self.sample % 50 == 0:
                force = self.actuator.step(command[self.sample // 50], stroke)
            else:
                force = self.actuator.force_at(stroke)
            self.sample += 1
            return force

        def force_at_stroke(self, stroke):
            return self.actuator.force_at(stroke)

    sea = jonswap_equal_energy_components(
        HYDRO, significant_height=.1042, peak_period=1.836,
        directions=np.array([0.]), spreading=np.array([1.]),
        gamma=3.3, phase=components[:, 3, None],
    )
    response = run_wavestar_published(
        HYDRO, sea, dt=.001, end_time=225, ramp_time=25,
        g=9.80665, pto_controller=PrescribedTorquePTO(), output_stride=50,
    )
    np.testing.assert_allclose(response.time, source_float[:, 0],
                               rtol=0, atol=1e-9)
    active = response.time >= 15
    source_force = source_pto[:, 33]
    peak_force = np.max(np.abs(source_force[active]))
    force_error = np.max(np.abs(response.pto_force - source_force))
    assert peak_force > 60
    assert force_error < .01 * peak_force
    errors = {
        "pitch_rad": np.max(np.abs(response.angle - source_float[:, 5])),
        "pitch_speed_rad_s": np.max(np.abs(
            response.angular_speed - source_float[:, 11])),
        "pto_stroke_m": np.max(np.abs(
            response.pto_stroke - source_pto[:, 3])),
        "axial_force_N": force_error,
    }
    assert np.isfinite(list(errors.values())).all()
    print("prescribed NMPC torque, independently advanced plant:", errors)

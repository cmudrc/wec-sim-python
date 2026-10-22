"""Keep the derived 225 s WaveStar split-step aggregate gate portable."""

from hashlib import sha256
from pathlib import Path

import numpy as np

from wecsim import (
    WaveStarNmpcActuator, WaveStarNmpcController, WaveStarNmpcObserver,
    WaveStarNmpcPredictor,
)
from wecsim.dynamics.wavestar import run_wavestar_published
from wecsim.environment.irregularWave import jonswap_equal_energy_components


ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests/fixtures/matlab/wavestar_nmpc_split_step_225s"
HYDRO = FIXTURE / "wavestar.h5"
PUBLISHED = ROOT / "tests/fixtures/matlab/wavestar_nmpc_full_step/baseline.npz"
SPLIT_STEP = FIXTURE / "source.npz"
SOURCE_SHA256 = "be50ad71195b5ff2d20c97bb3270845354c5a703d3829bf4aa0e915196ba7866"
HYDRO_SHA256 = "b2541ff01c319987d9157d1be9fbc9ec39a7660d0572d8e316c428b8aa3359d1"


def _check_source_fixture():
    assert sha256(SPLIT_STEP.read_bytes()).hexdigest() == SOURCE_SHA256
    assert sha256(HYDRO.read_bytes()).hexdigest() == HYDRO_SHA256


class _DelayedSampledPTO:
    """Hold the derived Simulink controller's one-sample actuator delay."""

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
            self.pending_command = self.previous_command = command
            self.commands.append(command)
        else:
            force = self.actuator.force_at(stroke)
        self.sample += 1
        return force

    def force_at_stroke(self, stroke):
        return self.actuator.force_at(stroke)


def test_full_split_step_controller_on_saved_stroke_and_prior_request():
    """The Python control law pairs source inputs before plant feedback."""
    _check_source_fixture()
    with np.load(SPLIT_STEP, allow_pickle=False) as source:
        time = source["time"]
        stroke = source["stroke"]
        request = source["command"]
    assert time.shape == stroke.shape == request.shape == (4501,)
    np.testing.assert_allclose(time, np.arange(4501) * .05,
                               rtol=0, atol=1e-9)
    assert np.isfinite(stroke).all() and np.isfinite(request).all()

    observer = WaveStarNmpcObserver(dt=.05)
    predictor = WaveStarNmpcPredictor(dt=.05)
    controller = WaveStarNmpcController(dt=.05)
    replayed = np.empty_like(request)
    for index, position in enumerate(stroke):
        prior = request[index - 1] if index else 0.
        state = observer.step(float(position), float(prior))
        forecast = predictor.step(float(state[-1]))
        replayed[index] = controller.step(state, forecast, float(time[index]))

    gap = replayed - request
    assert np.max(np.abs(gap)) < .01  # Published controller replay peak gate.
    assert np.sqrt(np.mean(gap**2)) < .0003  # Existing RMS gate.


def test_current_physical_path_matches_split_step_aggregates():
    """An aggregate pass does not promote the still-unpaired pointwise case."""
    _check_source_fixture()
    with np.load(PUBLISHED, allow_pickle=False) as published, \
         np.load(SPLIT_STEP, allow_pickle=False) as source:
        components = source["components"]
        assert components.shape == (500, 4)
        np.testing.assert_array_equal(components, published["components"])
        time = source["time"]
        np.testing.assert_allclose(time, np.arange(4501) * .05,
                                   rtol=0, atol=1e-9)
        assert all(source[key].shape == time.shape for key in (
            "pitch", "pitch_speed", "stroke", "pto_speed",
            "axial_force", "command",
        ))
        sea = jonswap_equal_energy_components(
            HYDRO, significant_height=.1042, peak_period=1.836,
            directions=np.array([0.]), spreading=np.array([1.]),
            gamma=3.3, phase=components[:, 3, None],
        )
        controller = _DelayedSampledPTO()
        response = run_wavestar_published(
            HYDRO, sea, dt=.001, end_time=225., ramp_time=25,
            g=9.80665, pto_controller=controller, output_stride=50,
        )
        np.testing.assert_allclose(response.time, time, rtol=0, atol=1e-9)
        command = np.asarray(controller.commands)
        assert command.shape == time.shape

        def rms(values):
            return float(np.sqrt(np.mean(np.square(values))))

        actual = {
            "pitch_rms": rms(response.angle),
            "pitch_peak": float(np.max(np.abs(response.angle))),
            "pitch_speed_rms": rms(response.angular_speed),
            "axial_force_rms": rms(response.pto_force),
            "axial_force_peak": float(np.max(np.abs(response.pto_force))),
            "requested_torque_rms": rms(command),
            "pto_work_J": float(np.trapezoid(
                response.pto_force * response.pto_speed, time,
            )),
        }
        expected = {
            "pitch_rms": rms(source["pitch"]),
            "pitch_peak": float(np.max(np.abs(source["pitch"]))),
            "pitch_speed_rms": rms(source["pitch_speed"]),
            "axial_force_rms": rms(source["axial_force"]),
            "axial_force_peak": float(np.max(np.abs(source["axial_force"]))),
            "requested_torque_rms": rms(source["command"]),
            "pto_work_J": float(np.trapezoid(
                source["axial_force"] * source["pto_speed"], time,
            )),
        }
        for label, value in expected.items():
            assert np.isfinite(value) and abs(value) > 1e-10, label
            gap = abs(actual[label] - value) / abs(value)
            assert np.isfinite(gap) and gap < .005, (
                f"{label}: {gap:.6%} exceeds the unchanged 0.5% limit"
            )

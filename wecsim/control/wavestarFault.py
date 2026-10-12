"""Sampled PTO controller from the published WaveStar fault application."""

import numpy as np

from wecsim.dynamics.friction import StribeckFriction
from wecsim.dynamics.wavestar import run_wavestar_published


class WaveStarFaultController:
    """Advance the published rotary damping controller one sample.

    ``noise`` and ``dropout`` describe a position-sensor disturbance. The
    actuator's second-order transfer function and the linkage conversion are
    those in the pinned fault application. ``step`` returns the axial PTO
    force in N; positive stroke shortens the B-to-C distance.
    """

    def __init__(self, *, gain=10., filter_frequency=30., sample_dt=.001):
        values = (gain, filter_frequency, sample_dt)
        if (not np.isfinite(values).all() or gain < 0
                or filter_frequency <= 0 or sample_dt <= 0):
            raise ValueError("WaveStar controller settings are invalid")
        self.gain = float(gain)
        self.filter_frequency = float(filter_frequency)
        self.sample_dt = float(sample_dt)
        self._previous_measurement = 0.
        self._previous_high = 0.
        self._previous_low = 0.
        self._previous_command = [0., 0.]
        self._previous_actuator = [0., 0.]

    @staticmethod
    def true_angle(stroke):
        length = .381408 + np.asarray(stroke, dtype=float)
        cosine = (length**2 - .412**2 - .2**2) / (-2 * .412 * .2)
        if not np.isfinite(cosine).all() or np.any(np.abs(cosine) > 1):
            raise ValueError("WaveStar PTO stroke is outside the linkage")
        return np.arccos(cosine) - 1.170165

    def step(self, stroke, *, noise=0., dropout=False):
        """Return axial PTO force for the current measured position sample."""
        if not np.isfinite([stroke, noise]).all():
            raise ValueError("WaveStar sensor input is invalid")
        true_angle = self.true_angle(stroke)
        measured = 0. if dropout else true_angle + noise
        angular_frequency = 2 * np.pi * self.filter_frequency
        rate = angular_frequency * self.sample_dt
        high = (2 * angular_frequency
                * (measured - self._previous_measurement)
                - (rate - 2) * self._previous_high) / (2 + rate)
        low = (rate * (high + self._previous_high)
               - (rate - 2) * self._previous_low) / (2 + rate)
        command = self.gain * low
        actuator = (
            .96294775 * command
            - 1.92342278 * self._previous_command[0]
            + .96053421 * self._previous_command[1]
            + 1.997777 * self._previous_actuator[0]
            - .99783206 * self._previous_actuator[1]
        )
        self._previous_measurement = measured
        self._previous_high = high
        self._previous_low = low
        self._previous_command = [command, self._previous_command[0]]
        self._previous_actuator = [actuator, self._previous_actuator[0]]
        angle = true_angle + 1.170165
        moment_arm = .412 * .2 * np.sin(angle) / (.381408 + stroke)
        return -actuator / moment_arm


def run_wavestar_fault_published(
    hydro_file, components, sensor_noise, sensor_dropout, *,
    end_time=141.2, output_stride=10, gain=10., filter_frequency=30.,
):
    """Integrate the pinned fault case using an exogenous sensor disturbance.

    Noise and dropout records contain the stochastic disturbance alone. The
    controller computes its true position from Python's evolving PTO stroke.
    MATLAB motion and PTO force are never inputs to the integration.
    """
    return run_wavestar_published(
        hydro_file, components, dt=.001, end_time=end_time,
        output_stride=output_stride,
        pto_controller=WaveStarFaultController(
            gain=gain, filter_frequency=filter_frequency,
        ),
        sensor_noise=sensor_noise, sensor_dropout=sensor_dropout,
        fault_joint_friction=StribeckFriction(.25, .1, .2, .001),
    )

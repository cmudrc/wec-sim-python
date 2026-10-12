"""Sampled three-coordinate PI force used by the WaveBot controller."""

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

from wecsim.hydrodynamics.bodyClass import BodyClass
from wecsim.dynamics.generalDynamics import BodyMotion, DynamicBody, GeneralizedDynamics
from wecsim.environment.irregularWave import IrregularComponents, synthesize_irregular_response


class WaveBotPIController:
    """Apply the published surge/heave/pitch gains to measured velocity.

    Gains are ``(Kp_surge, Ki_surge, Kp_heave, Ki_heave, Kp_pitch,
    Ki_pitch)``. The integrator advances after each output sample, matching
    the published Simulink discrete-time integrator.
    """

    def __init__(self, gains=(0., 0., 0., 0., 0., 0.)):
        self.integral = np.zeros(3)
        self.set_gains(gains)

    def set_gains(self, gains):
        values = np.asarray(gains, dtype=float)
        if values.shape != (6,) or not np.isfinite(values).all():
            raise ValueError("WaveBot PI gains need six finite values")
        self.gains = values.copy()

    def force(self, velocity):
        measured = np.asarray(velocity, dtype=float)
        if measured.shape != (3,) or not np.isfinite(measured).all():
            raise ValueError("WaveBot velocity needs three finite values")
        return self.gains[::2] * measured + self.gains[1::2] * self.integral

    def advance(self, velocity, dt):
        measured = np.asarray(velocity, dtype=float)
        if (measured.shape != (3,) or not np.isfinite(measured).all()
                or not np.isfinite(dt) or dt <= 0):
            raise ValueError("WaveBot integration needs velocity and positive dt")
        self.integral += dt * measured


class WaveBotPowerAutotuner:
    """Published power-only, fixed-weight WaveBot gain calculation."""

    def __init__(self, impedance: np.ndarray, frequency_hz: np.ndarray):
        impedance = np.asarray(impedance, dtype=complex)
        frequency = np.asarray(frequency_hz, dtype=float).ravel()
        if (impedance.shape != (3, 3, len(frequency))
                or len(frequency) < 2 or not np.isfinite(impedance).all()
                or not np.isfinite(frequency).all()
                or np.any(np.diff(frequency) <= 0)):
            raise ValueError("WaveBot impedance needs a finite 3x3 frequency series")
        self.omega = 2 * np.pi * np.arange(52, 257) / 256
        if (frequency[0] > self.omega[0] / (2 * np.pi)
                or frequency[-1] < self.omega[-1] / (2 * np.pi)):
            raise ValueError("WaveBot impedance does not cover 0.2 to 1 Hz")
        self.impedance = np.empty((205, 3, 3), dtype=complex)
        source_omega = 2 * np.pi * frequency
        for row in range(3):
            for col in range(3):
                entry = impedance[row, col]
                self.impedance[:, row, col] = (
                    np.interp(self.omega, source_omega, entry.real)
                    + 1j * np.interp(self.omega, source_omega, entry.imag)
                )
        self.surge_pitch_guess = np.array([-1000., 0., -200., 0.])

    def excitation_spectrum(self, velocity: np.ndarray,
                            delayed_command: np.ndarray) -> np.ndarray:
        """Estimate excitation from 0.25 s samples preceding an update."""
        velocity = np.asarray(velocity, dtype=float)
        delayed_command = np.asarray(delayed_command, dtype=float)
        if (velocity.ndim != 2 or velocity.shape[1:] != (3,)
                or delayed_command.shape != velocity.shape
                or not np.isfinite(velocity).all()
                or not np.isfinite(delayed_command).all()):
            raise ValueError("WaveBot FFT input needs matching finite three-channel samples")
        samples = min(len(velocity), 1024)
        data = np.zeros((1024, 6))
        if samples:
            data[-samples:, :3] = velocity[-samples:]
            data[-samples:, 3:] = delayed_command[-samples:]
        spectrum = np.fft.fft(np.hamming(1024)[:, None] * data,
                              axis=0)[52:257] / 1024
        return 2 * (
            np.einsum("fij,fj->fi", self.impedance, spectrum[:, :3])
            - spectrum[:, 3:]
        )

    def _heave_cost(self, gains, excitation):
        control = gains[0] - 1j * gains[1] / self.omega
        velocity = excitation / (self.impedance[:, 1, 1] - control)
        force = control * velocity
        kt, resistance, gear = 6.1745, .5, 12.4666
        term1 = (kt * 2 / 3 * gear * velocity
                 + resistance / (gear * kt) * force)
        term2 = force / (gear * kt)
        return 75 * float(np.real(np.vdot(term1, term2)))

    def _surge_pitch_cost(self, gains, excitation):
        control = np.column_stack((
            gains[0] - 1j * gains[1] / self.omega,
            gains[2] - 1j * gains[3] / self.omega,
        ))
        impedance = self.impedance[:, [0, 2]][:, :, [0, 2]].copy()
        impedance[:, 0, 0] -= control[:, 0]
        impedance[:, 1, 1] -= control[:, 1]
        try:
            velocity = np.linalg.solve(
                impedance, excitation[:, [0, 2], None],
            )[:, :, 0]
        except np.linalg.LinAlgError:
            return np.inf
        force = control * velocity
        kt, resistance = 6.1745, .5
        gear = np.array([12.4666, 3.])
        term1 = kt * 2 / 3 * gear * velocity + resistance / (gear * kt) * force
        term2 = force / (gear * kt)
        return 75 * float(np.real(np.sum(np.conj(term1) * term2)))

    def optimize(self, excitation: np.ndarray) -> np.ndarray:
        excitation = np.asarray(excitation, dtype=complex)
        if excitation.shape != (205, 3) or not np.isfinite(excitation).all():
            raise ValueError("WaveBot optimization needs 205 complex three-channel bins")
        options = {"maxiter": 10_000, "maxfev": 10_000,
                   "xatol": 1e-4, "fatol": 1e-4}
        heave = minimize(
            self._heave_cost, [-1000., 500.], args=(excitation[:, 1],),
            method="Nelder-Mead", options=options,
        )
        surge_pitch = minimize(
            self._surge_pitch_cost, self.surge_pitch_guess,
            args=(excitation,), method="Nelder-Mead", options=options,
        )
        if (not heave.success or not surge_pitch.success
                or not np.isfinite([*heave.x, *surge_pitch.x]).all()):
            raise RuntimeError("WaveBot power-only gain optimizer did not converge")
        self.surge_pitch_guess = surge_pitch.x.copy()
        return np.array([
            surge_pitch.x[0], surge_pitch.x[1], heave.x[0], heave.x[1],
            surge_pitch.x[2], surge_pitch.x[3],
        ])


@dataclass(frozen=True)
class WaveBotControlResponse:
    time: np.ndarray
    wave_elevation: np.ndarray
    body_position: np.ndarray
    body_velocity: np.ndarray
    actuator_command: np.ndarray
    gains: np.ndarray
    excitation_estimates: np.ndarray | None = None


def run_wavebot_prescribed_gains(
    hydro_file: str | Path,
    components: IrregularComponents,
    gains: np.ndarray,
    *,
    dt: float = .01,
    end_time: float = 100.,
    ramp_time: float = 20.,
    rho: float = 1025.,
    g: float = 9.81,
    source_linear_damping: bool = False,
) -> WaveBotControlResponse:
    """Advance the WaveBot using a supplied gain history.

    Supplying MATLAB's logged gains isolates plant, wave, and PI force
    parity. It does not independently reproduce the published autotuner.
    """
    return _run_wavebot_control(
        hydro_file, components, gains=gains, dt=dt, end_time=end_time,
        ramp_time=ramp_time, rho=rho, g=g,
        source_linear_damping=source_linear_damping,
    )


def run_wavebot_power_control(
    hydro_file: str | Path,
    components: IrregularComponents,
    impedance: np.ndarray,
    frequency_hz: np.ndarray,
    *,
    source_linear_damping: bool = False,
) -> WaveBotControlResponse:
    """Independently run the published 100 s adaptive power controller."""
    tuner = WaveBotPowerAutotuner(impedance, frequency_hz)
    return _run_wavebot_control(
        hydro_file, components, tuner=tuner,
        source_linear_damping=source_linear_damping,
    )


def _run_wavebot_control(
    hydro_file, components, *, gains=None, tuner=None, dt=.01,
    end_time=100., ramp_time=20., rho=1025., g=9.81,
    source_linear_damping=False,
) -> WaveBotControlResponse:
    if not isinstance(source_linear_damping, bool):
        raise TypeError("source_linear_damping must be a boolean")
    if (not np.isfinite([dt, end_time, ramp_time, rho, g]).all()
            or dt <= 0 or end_time <= 0 or ramp_time < 0
            or rho <= 0 or g <= 0):
        raise ValueError("WaveBot simulation settings are invalid")
    steps = round(end_time / dt)
    if not np.isclose(steps * dt, end_time, rtol=0, atol=1e-10):
        raise ValueError("WaveBot duration needs an integer number of steps")
    if (gains is None) == (tuner is None):
        raise ValueError("supply either prescribed gains or a power tuner")
    if tuner is None:
        gains = np.asarray(gains, dtype=float)
        if gains.shape != (steps + 1, 6) or not np.isfinite(gains).all():
            raise ValueError("gains must have six finite columns on the simulation grid")
    else:
        gains = np.zeros((steps + 1, 6))
    incident = synthesize_irregular_response(
        hydro_file, components, dt=dt, end_time=end_time,
        ramp_time=ramp_time, rho=rho, g=g,
    )
    body = BodyClass(str(hydro_file))
    body.bodyNumber = body.bodyTotal = 1
    body.readH5file()
    if int(np.asarray(body.dof).item()) != 6:
        raise ValueError("WaveBot HDF5 must contain one six-DOF body")
    body.mass = 1156.5
    zero = np.zeros((6, 6))
    body.hydroStiffness = zero.copy()
    body.viscDrag = {
        "Drag": zero.copy(), "cd": np.zeros(6),
        "characteristicArea": np.zeros(6),
    }
    body.linearDamping = zero.copy()
    memory_time = np.arange(round(10 / dt) + 1) * dt
    body.hydroForcePre(
        [], [0], len(memory_time), memory_time, [], dt, rho, g,
        "noWaveCIC", np.vstack((incident.time,
                                 np.zeros_like(incident.time))),
        1, 1, 0, 0, 0,
    )
    hydro = body.hydroForce
    mapping = np.zeros((6, 3))
    mapping[0, 0] = mapping[2, 1] = mapping[4, 2] = 1
    linear = np.diag([1000., 0., 1000., 0., 100., 0.])
    if source_linear_damping:
        linear = np.zeros((6, 6))
        linear[[0, 2, 4], 0] = [1000., 1000., 100.]
    drag = .5 * rho * np.array([1.15, 1.15, 1., .5, .5, 0.]) * np.array(
        [2.9568, 2.9568, 5.4739, 5.4739, 5.4739, 0.]
    )

    class AppliedForce:
        def __init__(self):
            self.controller = WaveBotPIController()
            self.previous_command = np.zeros(3)
            self.fft_velocity = []
            self.fft_command = []
            self.estimates = []

        def __call__(self, at_time, _q, speed):
            sample = int(np.clip(round(at_time / dt), 0, steps))
            if tuner is None:
                self.controller.set_gains(gains[sample])
            body_speed = mapping @ speed
            command = self.controller.force(speed)
            return (incident.excitation_force[sample]
                    - linear @ body_speed
                    - drag * np.abs(body_speed) * body_speed
                    + mapping @ (command * np.array([1., 1., -1.])))

        def commit_state(self, at_time, _q, speed):
            sample = int(round(at_time / dt))
            command = self.controller.force(speed)
            if tuner is not None:
                gains[sample] = self.controller.gains
                if sample % 25 == 0:
                    if sample % 800 == 0:
                        estimate = tuner.excitation_spectrum(
                            np.asarray(self.fft_velocity).reshape(-1, 3),
                            np.asarray(self.fft_command).reshape(-1, 3),
                        )
                        self.estimates.append(estimate)
                        if sample and sample < steps:
                            self.controller.set_gains(tuner.optimize(estimate))
                    self.fft_velocity.append(speed.copy())
                    self.fft_command.append(self.previous_command.copy())
            if tuner is not None and sample == 0:
                self.controller.set_gains([-1000., 0., -1000., 500., -200., 0.])
            self.controller.advance(speed, dt)
            self.previous_command = command

    def motion(q, _v):
        return BodyMotion(mapping @ q, mapping, np.zeros(6))

    applied = AppliedForce()
    dynamic = DynamicBody(
        rigid_mass=np.diag([1156.5] * 3 + [84.] * 3),
        added_mass=(np.asarray(hydro["fAddedMass"]),),
        damping=(zero,),
        restoring=np.asarray(hydro["linearHydroRestCoef"]),
        static_force=np.array([
            0., 0., (rho * float(np.asarray(body.dispVol).item()) - 1156.5) * g,
            0., 0., 0.,
        ]),
        reference_position=np.r_[np.asarray(body.cg, dtype=float).ravel(),
                                 np.zeros(3)],
        motion=motion,
        excitation=lambda _time: np.zeros(6),
        radiation_kernel=np.asarray(hydro["irkb"]),
        state_excitation=applied,
    )
    system = GeneralizedDynamics(
        (dynamic,), 3,
        pto_stiffness=np.diag([24_000., 5_000., 1_000.]),
        pto_damping=np.diag([5., 5., 0.]),
    )
    response = system.integrate(dt=dt, end_time=end_time)
    controller = WaveBotPIController()
    actuator = np.empty((steps + 1, 3))
    for sample, speed in enumerate(response.speed):
        controller.set_gains(gains[sample])
        actuator[sample] = controller.force(speed)
        controller.advance(speed, dt)
    return WaveBotControlResponse(
        response.time, incident.elevation,
        response.body_position[:, 0], response.body_velocity[:, 0],
        actuator, gains,
        (np.asarray(applied.estimates) if tuner is not None else None),
    )

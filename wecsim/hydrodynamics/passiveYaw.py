"""Directional excitation for a body that yaws with respect to the waves."""

from dataclasses import dataclass

import numpy as np
from scipy.interpolate import CubicSpline, make_interp_spline

from wecsim.environment.irregularWave import IrregularComponents


@dataclass(frozen=True)
class PassiveYawExcitation:
    """Regular-wave force from BEM headings at the body's current yaw angle.

    The hydrodynamic force is interpolated at incident heading minus body yaw
    and rotated from body axes into world axes. Heading interpolation is
    continuous; the pinned MATLAB block instead reuses coefficients until yaw
    changes by its configured threshold.
    """

    directions: np.ndarray
    real: np.ndarray  # six DOFs by heading, in N/m or Nm/m
    imaginary: np.ndarray
    incident_direction: float
    omega: float
    amplitude: float
    ramp_time: float

    @classmethod
    def from_hydro_data(cls, hydro_data, *, omega, incident_direction,
                        amplitude, ramp_time, rho, g,
                        spline_frequency=False):
        parameters = hydro_data["simulation_parameters"]
        headings = np.asarray(parameters["wave_dir"], dtype=float).ravel()
        frequency = np.asarray(parameters["w"], dtype=float).ravel()
        if (headings.size < 3 or frequency.size < 2
                or not np.isfinite(headings).all()
                or not np.isfinite(frequency).all()
                or not np.all(np.diff(headings) > 0)
                or not np.all(np.diff(frequency) > 0)
                or np.max(np.diff(np.r_[headings, headings[0] + 360])) >= 180
                or omega < frequency[0] or omega > frequency[-1]):
            raise ValueError(
                "passive yaw needs full-circle BEM headings and an in-range frequency"
            )
        excitation = hydro_data["hydro_coeffs"]["excitation"]
        real = np.asarray(excitation["re"], dtype=float)
        imaginary = np.asarray(excitation["im"], dtype=float)
        if (real.shape != (6, headings.size, frequency.size)
                or imaginary.shape != real.shape
                or not np.isfinite(real).all() or not np.isfinite(imaginary).all()):
            raise ValueError("passive-yaw excitation coefficients are incomplete")
        if spline_frequency:
            # MATLAB's single-direction bank calls interp1(..., 'spline');
            # its full-direction path uses linear interp2 instead.
            real = CubicSpline(frequency, real, axis=2)(omega) * rho * g
            imaginary = CubicSpline(frequency, imaginary, axis=2)(omega) * rho * g
        else:
            real = make_interp_spline(frequency, real, k=1, axis=2)(omega) * rho * g
            imaginary = make_interp_spline(
                frequency, imaginary, k=1, axis=2,
            )(omega) * rho * g
        return cls(headings, real, imaginary, incident_direction,
                   omega, amplitude, ramp_time)

    def force(self, time: float, yaw: float, *,
              coefficient_heading: float | None = None) -> np.ndarray:
        """Return the six-component excitation in the world frame."""
        relative_heading = (
            self.incident_direction - np.degrees(yaw)
            if coefficient_heading is None else coefficient_heading
        ) % 360
        real = np.array([
            np.interp(relative_heading, self.directions, row, period=360)
            for row in self.real
        ])
        imaginary = np.array([
            np.interp(relative_heading, self.directions, row, period=360)
            for row in self.imaginary
        ])
        ramp = (1.0 if self.ramp_time == 0 or time >= self.ramp_time
                else (1 - np.cos(np.pi * time / self.ramp_time)) / 2)
        local = self.amplitude * ramp * (
            real * np.cos(self.omega * time)
            - imaginary * np.sin(self.omega * time)
        )
        world = local.copy()
        c, s = np.cos(yaw), np.sin(yaw)
        for first in (0, 3):
            world[first] = c * local[first] - s * local[first + 1]
            world[first + 1] = s * local[first] + c * local[first + 1]
        return world


@dataclass(frozen=True)
class NearestHeadingExcitation:
    """Select the nearest BEM heading as in the variable-hydro yaw example.

    The published direction-bank files share mass, restoring, and radiation
    coefficients; their excitation coefficients differ by heading. Unlike
    passive yaw, the variable-hydro source does not rotate that selected
    excitation by body yaw. This selector uses the full-direction HDF5 input.
    """

    model: PassiveYawExcitation
    headings: np.ndarray

    def __post_init__(self):
        headings = np.asarray(self.headings, dtype=float)
        if (headings.ndim != 1 or headings.size < 2
                or not np.isfinite(headings).all()
                or not np.all(np.diff(headings) > 0)
                or headings[0] < -180 or headings[-1] > 180):
            raise ValueError("heading bank needs ordered directions in [-180, 180]")
        object.__setattr__(self, "headings", headings)

    def heading(self, yaw: float) -> float:
        relative = self.model.incident_direction - np.degrees(yaw)
        return float(self.headings[np.argmin(np.abs(self.headings - relative))])

    def force(self, time: float, yaw: float) -> np.ndarray:
        return self.model.force(
            time, 0.0, coefficient_heading=self.heading(yaw),
        )


@dataclass(frozen=True)
class SampledPassiveYawExcitation:
    """Broadband force at each BEM heading, interpolated at the body's yaw.

    Frequency interpolation and the irregular realization are prepared once.
    The dynamics can then evaluate changing yaw during a sampled radiation
    step without generating new random phases or recomputing the spectrum.
    """

    time: np.ndarray
    directions: np.ndarray
    incident_directions: np.ndarray
    force_grid: np.ndarray  # time, incident direction, BEM heading, six DOFs
    elevation: np.ndarray

    @classmethod
    def from_hydro_data(cls, hydro_data, components: IrregularComponents,
                        *, dt: float, end_time: float, ramp_time: float,
                        rho: float, g: float):
        if not isinstance(components, IrregularComponents):
            raise TypeError("components must be an irregular-wave realization")
        if (not np.isfinite([dt, end_time, ramp_time, rho, g]).all()
                or dt <= 0 or end_time < 0 or ramp_time < 0
                or rho <= 0 or g <= 0):
            raise ValueError("time and fluid parameters must be finite and valid")
        count = round(end_time / dt) + 1
        if not np.isclose((count - 1) * dt, end_time, rtol=0, atol=1e-10):
            raise ValueError("end_time must be a multiple of dt")
        omega = np.asarray(components.omega, dtype=float).ravel()
        amplitude = np.asarray(components.spectral_amplitude, dtype=float).ravel()
        width = np.asarray(components.d_omega, dtype=float).ravel()
        incident = np.asarray(components.directions, dtype=float).ravel()
        spread = np.asarray(components.spreading, dtype=float).ravel()
        phase = np.asarray(components.phase, dtype=float)
        if (len(omega) < 2 or len(incident) < 1
                or amplitude.shape != omega.shape or width.shape != omega.shape
                or spread.shape != incident.shape
                or phase.shape != (len(omega), len(incident))
                or not all(np.isfinite(values).all() for values in
                           (omega, amplitude, width, incident, spread, phase))
                or np.any(np.diff(omega) <= 0) or np.any(amplitude < 0)
                or np.any(width <= 0) or np.any(spread < 0)
                or not np.isclose(spread.sum(), 1, rtol=0, atol=1e-10)):
            raise ValueError("irregular-wave components are inconsistent")
        parameters = hydro_data["simulation_parameters"]
        headings = np.asarray(parameters["wave_dir"], dtype=float).ravel()
        frequency = np.asarray(parameters["w"], dtype=float).ravel()
        if (headings.size < 3 or frequency.size < 2
                or not np.isfinite(headings).all()
                or not np.isfinite(frequency).all()
                or not np.all(np.diff(headings) > 0)
                or not np.all(np.diff(frequency) > 0)
                or np.max(np.diff(np.r_[headings, headings[0] + 360])) >= 180
                or omega[0] < frequency[0] or omega[-1] > frequency[-1]):
            raise ValueError("irregular passive yaw needs full-circle BEM headings and frequencies")
        excitation = hydro_data["hydro_coeffs"]["excitation"]

        def frequency_samples(field):
            raw = np.asarray(field, dtype=float)
            if (raw.shape != (6, len(headings), len(frequency))
                    or not np.isfinite(raw).all()):
                raise ValueError("irregular passive-yaw excitation is incomplete")
            return (CubicSpline(frequency, raw, axis=2)(omega)
                    .transpose(2, 1, 0).reshape(len(omega), -1) * rho * g)

        real = frequency_samples(excitation["re"])
        imaginary = frequency_samples(excitation["im"])
        time = np.arange(count) * dt
        forces = np.empty((count, len(incident), len(headings), 6))
        elevation = np.zeros(count)
        for sea in range(len(incident)):
            height = np.sqrt(amplitude * width * spread[sea])
            for start in range(0, count, 128):
                stop = min(start + 128, count)
                angle = (time[start:stop, None] * omega[None, :]
                         + phase[None, :, sea])
                cosine = np.cos(angle) * height[None, :]
                sine = np.sin(angle) * height[None, :]
                forces[start:stop, sea] = (
                    cosine @ real - sine @ imaginary
                ).reshape(stop - start, len(headings), 6)
                elevation[start:stop] += cosine.sum(axis=1)
        ramp = np.ones(count)
        if ramp_time > 0:
            early = time < ramp_time
            ramp[early] = (1 - np.cos(np.pi * time[early] / ramp_time)) / 2
        return cls(time, headings, incident, forces * ramp[:, None, None, None],
                   elevation * ramp)

    def force(self, at_time: float, yaw: float, *, coefficient_yaw: float | None = None) -> np.ndarray:
        """Return world-frame excitation at one sample and yaw angle."""
        dt = self.time[1] - self.time[0] if len(self.time) > 1 else None
        index = round(at_time / dt) if dt is not None else 0
        if (index < 0 or index >= len(self.time)
                or not np.isclose(self.time[index], at_time, rtol=0, atol=1e-8)):
            raise ValueError("sampled passive-yaw force needs a grid time")
        local = np.zeros(6)
        if coefficient_yaw is None:
            coefficient_yaw = yaw
        start = self.directions[0]
        for sea, heading in enumerate(self.incident_directions):
            relative = ((heading - np.degrees(coefficient_yaw) - start) % 360) + start
            left = np.searchsorted(self.directions, relative, side="right") - 1
            right = (left + 1) % len(self.directions)
            next_heading = (self.directions[right] if right > left
                            else self.directions[0] + 360)
            fraction = (relative - self.directions[left]) / (
                next_heading - self.directions[left]
            )
            local += ((1 - fraction) * self.force_grid[index, sea, left]
                      + fraction * self.force_grid[index, sea, right])
        world = local.copy()
        c, s = np.cos(yaw), np.sin(yaw)
        for first in (0, 3):
            world[first] = c * local[first] - s * local[first + 1]
            world[first + 1] = s * local[first] + c * local[first + 1]
        return world


@dataclass(frozen=True)
class NearestSampledHeadingExcitation:
    """Select an irregular-wave force bank without passive-yaw rotation."""

    model: SampledPassiveYawExcitation
    headings: np.ndarray

    def __post_init__(self):
        headings = np.asarray(self.headings, dtype=float)
        if (len(self.model.incident_directions) != 1
                or headings.ndim != 1 or headings.size < 2
                or not np.isfinite(headings).all()
                or not np.all(np.diff(headings) > 0)
                or headings[0] < -180 or headings[-1] > 180):
            raise ValueError("sampled heading bank needs one sea and ordered directions in [-180, 180]")
        object.__setattr__(self, "headings", headings)

    def heading(self, yaw: float) -> float:
        relative = self.model.incident_directions[0] - np.degrees(yaw)
        return float(self.headings[np.argmin(np.abs(self.headings - relative))])

    def force(self, time: float, yaw: float) -> np.ndarray:
        incident = self.model.incident_directions[0]
        coefficient_yaw = np.deg2rad(incident - self.heading(yaw))
        return self.model.force(time, 0.0, coefficient_yaw=coefficient_yaw)

    def __call__(self, at_time: float, coordinate: np.ndarray,
                 speed: np.ndarray) -> np.ndarray:
        return self.force(at_time, float(coordinate[0]))

    def freeze_step_excitation(
        self, at_time: float, accepted_coordinate: np.ndarray,
        trial_coordinate: np.ndarray, next_coordinate: np.ndarray,
    ) -> np.ndarray | None:
        """Use the accepted heading only when endpoint trials cross a bank edge."""
        if self.heading(float(trial_coordinate[0])) == self.heading(float(next_coordinate[0])):
            return None
        return self.force(at_time, float(accepted_coordinate[0]))


class HeldPassiveYawExcitation:
    """Apply the source's sampled heading-coefficient update threshold.

    Force evaluation during a trial step is side-effect free. The integrator
    commits a heading only after it accepts the corresponding body state.
    """

    def __init__(self, model: SampledPassiveYawExcitation, threshold: float):
        if (len(model.incident_directions) != 1 or not np.isfinite(threshold)
                or threshold <= 0):
            raise ValueError("held passive yaw needs one incident direction and positive threshold")
        self.model = model
        self.threshold = threshold
        self.last_heading: float | None = None
        self.force_history: list[np.ndarray] = []

    def _heading(self, yaw: float) -> float:
        relative = float(self.model.incident_directions[0] - np.degrees(yaw))
        directions = self.model.directions
        wrapped_distance = (relative - directions + 180) % 360 - 180
        nearest_index = int(np.argmin(np.abs(wrapped_distance)))
        if abs(wrapped_distance[nearest_index]) <= self.threshold:
            # The pinned MATLAB block selects a tabulated BEM heading when
            # the current relative heading lies inside its threshold.
            return relative - wrapped_distance[nearest_index]
        if (self.last_heading is None
                or abs(relative - self.last_heading) > self.threshold):
            return relative
        return self.last_heading

    def _force(self, at_time: float, yaw: float, heading: float) -> np.ndarray:
        coefficient_yaw = np.deg2rad(self.model.incident_directions[0] - heading)
        return self.model.force(at_time, yaw, coefficient_yaw=coefficient_yaw)

    def __call__(self, at_time: float, coordinate: np.ndarray,
                 speed: np.ndarray) -> np.ndarray:
        yaw = float(coordinate[0])
        return self._force(at_time, yaw, self._heading(yaw))

    def commit(self, at_time: float, coordinate: np.ndarray) -> None:
        yaw = float(coordinate[0])
        heading = self._heading(yaw)
        self.force_history.append(self._force(at_time, yaw, heading))
        self.last_heading = heading

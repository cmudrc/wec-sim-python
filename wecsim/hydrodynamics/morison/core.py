"""Morison element geometry and shared wave/current kinematics."""

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class MorisonElement:
    point: tuple[float, float, float]  # body-local, relative to center of gravity
    drag_coefficient: tuple[float, float, float]
    added_mass_coefficient: tuple[float, float, float]
    area: tuple[float, float, float]
    volume: float
    phase_mode: str = "directional"  # "matlab_shared" replays the pinned source loop
    axis: tuple[float, float, float] | None = None  # body-local option-2 axis


@dataclass(frozen=True)
class FixedMorisonResponse:
    time: np.ndarray
    wave_elevation: np.ndarray
    force: np.ndarray  # physical force and moment at the body's center of gravity


def finite_depth_wavenumber(
    omega: np.ndarray, *, water_depth: float, gravity: float = 9.81,
) -> np.ndarray:
    """Use WEC-Sim's 100 fixed-point dispersion iterations."""
    frequency = np.asarray(omega, dtype=float)
    if (frequency.ndim != 1 or len(frequency) == 0
            or not np.isfinite(frequency).all() or np.any(frequency <= 0)
            or not np.isfinite([water_depth, gravity]).all()
            or water_depth <= 0 or gravity <= 0):
        raise ValueError("wavenumber needs positive frequency, depth, and gravity")
    k = frequency**2 / gravity
    for _ in range(100):
        k = frequency**2 / gravity / np.tanh(k * water_depth)
    return k


def _validate_current(speed, direction, profile, depth):
    if (not np.isfinite([speed, direction]).all() or speed < 0
            or not -360 <= direction <= 360
            or profile not in ("uniform", "power", "linear")
            or (profile != "uniform" and
                (depth is None or not np.isfinite(depth) or depth <= 0))):
        raise ValueError("regular Morison current settings are invalid")


def _current_speed_at_depth(z, time, ramp_time, speed, profile, depth):
    if not speed:
        return 0.0
    if profile != "uniform":
        if z <= -depth:
            return 0.0
        fraction = 1 + z / depth
        speed *= fraction ** (1 / 7 if profile == "power" else 1)
    if ramp_time and time < ramp_time:
        speed *= (1 - np.cos(np.pi * time / ramp_time)) / 2
    return speed

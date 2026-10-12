"""Wave and current forces on stationary Morison elements."""

from typing import Sequence

import numpy as np

from wecsim.environment.irregularWave import IrregularComponents

from .core import MorisonElement, FixedMorisonResponse, finite_depth_wavenumber


def solve_fixed_morison_irregular(
    components: IrregularComponents,
    elements: Sequence[MorisonElement],
    *,
    center_gravity: Sequence[float],
    water_depth: float,
    dt: float,
    end_time: float,
    ramp_time: float,
    rho: float = 1025.0,
    g: float = 9.81,
    current_speed: float = 0.0,
    current_direction: float = 0.0,
    current_profile: str = "uniform",
    current_depth: float | None = None,
) -> FixedMorisonResponse:
    """Evaluate a stationary body's Cartesian Morison force at every sample."""
    center = np.asarray(center_gravity, dtype=float)
    if center.shape != (3,) or not np.isfinite(center).all():
        raise ValueError("center_gravity needs three finite coordinates")
    if (not np.isfinite([dt, end_time, ramp_time, rho, g]).all()
            or dt <= 0 or end_time < 0 or ramp_time < 0 or rho <= 0 or g <= 0):
        raise ValueError("fixed Morison time and fluid settings are invalid")
    if (not np.isfinite([current_speed, current_direction]).all()
            or current_speed < 0 or not -360 <= current_direction <= 360
            or current_profile not in ("uniform", "power", "linear")
            or (current_profile != "uniform" and
                (current_depth is None or not np.isfinite(current_depth)
                 or current_depth <= 0))):
        raise ValueError("fixed Morison current settings are invalid")
    steps = round(end_time / dt)
    if not np.isclose(steps * dt, end_time, rtol=0, atol=1e-10):
        raise ValueError("end_time must be a multiple of dt")
    if not elements:
        raise ValueError("a fixed Morison body needs at least one element")
    prepared = []
    for element in elements:
        if not isinstance(element, MorisonElement):
            raise TypeError("elements must be MorisonElement values")
        point = np.asarray(element.point, dtype=float)
        drag = np.asarray(element.drag_coefficient, dtype=float)
        added = np.asarray(element.added_mass_coefficient, dtype=float)
        area = np.asarray(element.area, dtype=float)
        if (any(value.shape != (3,) or not np.isfinite(value).all()
                for value in (point, drag, added, area))
                or np.any(drag < 0) or np.any(added < 0) or np.any(area < 0)
                or not np.isfinite(element.volume) or element.volume <= 0
                or element.phase_mode not in ("directional", "matlab_shared")
                or element.axis is not None):
            raise ValueError("Morison point, coefficients, area, and volume are invalid")
        prepared.append((point, drag, added, area, float(element.volume),
                         element.phase_mode))

    omega = np.asarray(components.omega, dtype=float)
    spectrum = np.asarray(components.spectral_amplitude, dtype=float)
    width = np.asarray(components.d_omega, dtype=float)
    directions = np.asarray(components.directions, dtype=float)
    spreading = np.asarray(components.spreading, dtype=float)
    phase = np.asarray(components.phase, dtype=float)
    if (omega.ndim != 1 or len(omega) < 2
            or spectrum.shape != omega.shape or width.shape != omega.shape
            or directions.ndim != 1 or spreading.shape != directions.shape
            or phase.shape != (len(omega), len(directions))
            or any(not np.isfinite(value).all() for value in
                   (omega, spectrum, width, directions, spreading, phase))
            or np.any(omega <= 0) or np.any(spectrum < 0)
            or np.any(width <= 0) or np.any(spreading < 0)
            or not np.isclose(spreading.sum(), 1, atol=1e-12)):
        raise ValueError("invalid directional irregular-wave components")
    if current_speed and len(directions) != 1:
        raise ValueError("fixed Morison current currently needs one incident heading")
    k = finite_depth_wavenumber(omega, water_depth=water_depth, gravity=g)
    time = np.arange(steps + 1, dtype=float) * dt
    ramp = np.ones_like(time)
    if ramp_time:
        before = time < ramp_time
        ramp[before] = (1 - np.cos(np.pi * time[before] / ramp_time)) / 2
    elevation = np.zeros_like(time)
    force = np.zeros((len(time), 6))
    headings = np.deg2rad(directions)
    current_heading = np.deg2rad(current_direction)
    current_axis = np.array([np.cos(current_heading), np.sin(current_heading), 0.0])

    for start in range(0, len(time), 512):
        stop = min(start + 512, len(time))
        t = time[start:stop, None]
        local_ramp = ramp[start:stop]
        for heading_index, heading in enumerate(headings):
            direction = np.array([np.cos(heading), np.sin(heading)])
            amplitude = np.sqrt(spreading[heading_index] * spectrum * width)
            origin_phase = t * omega + phase[:, heading_index]
            elevation[start:stop] += local_ramp * (
                np.cos(origin_phase) @ amplitude
            )
            for point, drag, added, area, volume, phase_mode in prepared:
                world = center + point
                if world[2] > 0:
                    continue
                kh = k * water_depth
                kz = k * world[2]
                horizontal = np.where(
                    kh > np.pi, np.exp(kz),
                    np.cosh(kz + kh) / np.cosh(kh),
                )
                vertical = np.where(
                    kh > np.pi, np.exp(kz),
                    np.sinh(kz + kh) / np.cosh(kh),
                )
                # The pinned irregWaveMorison.m indexes randPhase(jj,1)
                # inside its heading loop. Keep that source behavior opt-in.
                force_phase = (phase[:, 0] if phase_mode == "matlab_shared"
                               else phase[:, heading_index])
                argument = (t * omega + force_phase
                            - k * (world[:2] @ direction))
                cosine = np.cos(argument)
                sine = np.sin(argument)
                horizontal_velocity = local_ramp * (
                    cosine @ (amplitude * horizontal * g * k / omega)
                )
                vertical_velocity = -local_ramp * (
                    sine @ (amplitude * vertical * g * k / omega)
                )
                horizontal_acceleration = -local_ramp * (
                    sine @ (amplitude * horizontal * g * k)
                )
                vertical_acceleration = -local_ramp * (
                    cosine @ (amplitude * vertical * g * k)
                )
                velocity = np.column_stack((
                    horizontal_velocity * direction[0],
                    horizontal_velocity * direction[1],
                    vertical_velocity,
                ))
                if current_speed:
                    if current_profile == "uniform":
                        current_at_point = current_speed
                    elif world[2] <= -current_depth:
                        current_at_point = 0.0
                    else:
                        fraction = 1 + world[2] / current_depth
                        exponent = 1 / 7 if current_profile == "power" else 1
                        current_at_point = current_speed * fraction**exponent
                    velocity += local_ramp[:, None] * current_at_point * current_axis
                acceleration = np.column_stack((
                    horizontal_acceleration * direction[0],
                    horizontal_acceleration * direction[1],
                    vertical_acceleration,
                ))
                contribution = (rho * volume * (1 + added) * acceleration
                                + 0.5 * rho * drag * area * velocity
                                * np.abs(velocity))
                force[start:stop, :3] += contribution
                force[start:stop, 3:] += np.cross(point, contribution)
    return FixedMorisonResponse(time, elevation, force)

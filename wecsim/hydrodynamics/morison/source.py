"""Source-specific Morison force reconstructions used in paired audits."""

from typing import Sequence

import numpy as np

from wecsim.environment.irregularWave import IrregularComponents

from .core import (
    MorisonElement, finite_depth_wavenumber,
    _current_speed_at_depth, _validate_current,
)


def regular_morison_source_force(
    elements: Sequence[MorisonElement], *, time: float,
    position: Sequence[float], velocity: Sequence[float],
    acceleration: Sequence[float], wave_height: float, wave_period: float,
    direction: float, water_depth: float, ramp_time: float,
    rho: float = 1025.0, g: float = 9.81,
    current_speed: float = 0.0, current_direction: float = 0.0,
    current_profile: str = "uniform", current_depth: float | None = None,
    body_morison: int = 1,
    element_axes: Sequence[Sequence[float]] | None = None,
) -> np.ndarray:
    """Evaluate pinned ``regWaveMorison.m`` at a moving-body state.

    The six state components are surge, sway, heave, roll, pitch, and yaw.
    This source-law diagnostic does not solve a coupled WEC trajectory. It
    retains the source's rotation and local-point
    angular kinematics so that a comparison can expose source-specific
    behavior. Option 2 requires one body-local axis per element.
    """
    state = [np.asarray(value, dtype=float) for value in
             (position, velocity, acceleration)]
    if any(value.shape != (6,) or not np.isfinite(value).all()
           for value in state):
        raise ValueError("Morison position, velocity, and acceleration need finite six-vectors")
    if (not np.isfinite([time, wave_height, wave_period, direction,
                         water_depth, ramp_time, rho, g]).all()
            or time < 0 or wave_height < 0 or wave_period <= 0
            or water_depth <= 0 or ramp_time < 0 or rho <= 0 or g <= 0):
        raise ValueError("regular Morison wave and fluid settings are invalid")
    _validate_current(current_speed, current_direction,
                      current_profile, current_depth)
    if not elements:
        raise ValueError("regular Morison force needs at least one element")
    if body_morison not in (1, 2):
        raise ValueError("body_morison must be Cartesian (1) or normal/tangential (2)")
    if body_morison == 1 and element_axes is not None:
        raise ValueError("element_axes are only used by normal/tangential mode")
    if body_morison == 2:
        if element_axes is None or len(element_axes) != len(elements):
            raise ValueError("normal/tangential mode needs one axis per element")
        for axis in element_axes:
            vector = np.asarray(axis, dtype=float)
            if (vector.shape != (3,) or not np.isfinite(vector).all()
                    or np.linalg.norm(vector) == 0):
                raise ValueError("Morison element axes need finite nonzero three-vectors")

    pose, speed, accel = state
    roll, pitch, yaw = pose[3:]
    c4, s4 = np.cos(roll), np.sin(roll)
    c5, s5 = np.cos(pitch), np.sin(pitch)
    c6, s6 = np.cos(yaw), np.sin(yaw)
    # Keep the pinned MATLAB source's matrix, including its first-row term.
    rotation = np.array([
        [c5 * c6, c4 * s6 + s4 * s5 * c6, s4 * s6 - c4 * s5 * s6],
        [-c5 * s6, c4 * c6 - s4 * s5 * s6, s4 * c6 + c4 * s5 * s6],
        [s5, -s4 * c5, c4 * c5],
    ])
    omega = 2 * np.pi / wave_period
    k = finite_depth_wavenumber(
        np.array([omega]), water_depth=water_depth, gravity=g,
    )[0]
    heading = np.deg2rad(direction)
    wave_axis = np.array([np.cos(heading), np.sin(heading)])
    amplitude = wave_height / 2
    if ramp_time and time < ramp_time:
        amplitude *= (1 - np.cos(np.pi * time / ramp_time)) / 2
    result = np.zeros(6)
    for index, element in enumerate(elements):
        if not isinstance(element, MorisonElement):
            raise TypeError("elements must be MorisonElement values")
        point, cd, ca, area = (
            np.asarray(value, dtype=float) for value in (
                element.point, element.drag_coefficient,
                element.added_mass_coefficient, element.area,
            )
        )
        if (any(value.shape != (3,) or not np.isfinite(value).all()
                for value in (point, cd, ca, area))
                or np.any(cd < 0) or np.any(ca < 0) or np.any(area < 0)
                or not np.isfinite(element.volume) or element.volume <= 0):
            raise ValueError("Morison element coefficients and geometry are invalid")
        rotated_point = rotation @ point
        world = pose[:3] + rotated_point
        if world[2] > 0:
            continue
        angular_cross_point = np.cross(speed[3:], point)
        body_velocity = speed[:3] + angular_cross_point
        body_acceleration = (accel[:3] + np.cross(accel[3:], point)
                             + np.cross(speed[3:], angular_cross_point))
        kh, kz = k * water_depth, k * world[2]
        if kh > np.pi:
            horizontal = vertical = np.exp(kz)
        else:
            horizontal = np.cosh(kz + kh) / np.cosh(kh)
            vertical = np.sinh(kz + kh) / np.cosh(kh)
        phase = omega * time - k * (world[:2] @ wave_axis)
        horizontal_velocity = amplitude * horizontal * np.cos(phase) * g * k / omega
        vertical_velocity = -amplitude * vertical * np.sin(phase) * g * k / omega
        horizontal_acceleration = -amplitude * horizontal * np.sin(phase) * g * k
        vertical_acceleration = -amplitude * vertical * np.cos(phase) * g * k
        fluid_velocity = np.r_[horizontal_velocity * wave_axis, vertical_velocity]
        fluid_velocity[:2] += _current_speed_at_depth(
            world[2], time, ramp_time, current_speed, current_profile,
            current_depth,
        ) * np.array([np.cos(np.deg2rad(current_direction)),
                      np.sin(np.deg2rad(current_direction))])
        fluid_acceleration = np.r_[horizontal_acceleration * wave_axis,
                                   vertical_acceleration]
        if body_morison == 2:
            axis = rotation @ np.asarray(element_axes[index], dtype=float)
            axis_norm = np.linalg.norm(axis)
            if axis_norm == 0:
                raise ValueError("rotated Morison axis has zero length")
            axis_norm_squared = axis_norm**2
            fluid_velocity_t = axis * (axis @ fluid_velocity / axis_norm_squared)
            body_velocity_t = axis * (axis @ body_velocity / axis_norm_squared)
            fluid_acceleration_t = axis * (
                axis @ fluid_acceleration / axis_norm_squared
            )
            body_acceleration_t = axis * (
                axis @ body_acceleration / axis_norm_squared
            )
            relative_t = fluid_velocity_t - body_velocity_t
            relative_n = ((fluid_velocity - fluid_velocity_t)
                          - (body_velocity - body_velocity_t))
            force = (0.5 * rho * (
                cd[0] * area[0] * relative_n * np.linalg.norm(relative_n)
                + cd[1] * area[1] * relative_t * np.linalg.norm(relative_t)
            ) + rho * element.volume * (
                fluid_acceleration
                + ca[0] * ((fluid_acceleration - fluid_acceleration_t)
                           - (body_acceleration - body_acceleration_t))
                + ca[1] * (fluid_acceleration_t - body_acceleration_t)
            ))
        else:
            relative_velocity = fluid_velocity - body_velocity
            # These row-vector coefficient transforms are the source block's
            # option-1 convention; they are not a rotated drag tensor.
            area_rot = np.abs(area @ rotation)
            cd_rot = np.abs(cd) @ rotation
            ca_rot = np.abs(ca @ rotation)
            force = (0.5 * rho * cd_rot * area_rot * relative_velocity
                     * np.abs(relative_velocity)
                     + rho * element.volume * (
                         fluid_acceleration
                         + ca_rot * (fluid_acceleration - body_acceleration)))
        result[:3] += force
        result[3:] += np.cross(rotated_point, force)
    return result


def irregular_morison_source_drag(
    elements: Sequence[MorisonElement], *, time: float,
    position: Sequence[float], velocity: Sequence[float],
    components: IrregularComponents, water_depth: float, ramp_time: float,
    rho: float = 1000.0, g: float = 9.81,
) -> np.ndarray:
    """Replay the pinned single-heading Cartesian drag law at one body state.

    This is a source diagnostic for the published Desalination flap's
    zero-volume Morison elements. It retains ``irregWaveMorison.m``'s
    rotation and local-point velocity conventions; physical body motion
    need not use those conventions.
    """
    pose = np.asarray(position, dtype=float)
    speed = np.asarray(velocity, dtype=float)
    omega = np.asarray(components.omega, dtype=float).ravel()
    spectrum = np.asarray(components.spectral_amplitude, dtype=float).ravel()
    width = np.asarray(components.d_omega, dtype=float).ravel()
    phase = np.asarray(components.phase, dtype=float)
    heading = np.asarray(components.directions, dtype=float).ravel()
    spread = np.asarray(components.spreading, dtype=float).ravel()
    if (pose.shape != (6,) or speed.shape != (6,)
            or heading.shape != (1,) or spread.shape != (1,)
            or phase.shape != (len(omega), 1)
            or spectrum.shape != omega.shape or width.shape != omega.shape
            or not np.isfinite(pose).all() or not np.isfinite(speed).all()
            or not np.isfinite(phase).all()
            or not np.isfinite(spectrum).all()
            or not np.isfinite(width).all()
            or not np.isfinite(heading).all()
            or not np.isfinite(spread).all()
            or np.any(spectrum < 0) or np.any(width <= 0)
            or not np.isfinite([time, water_depth, ramp_time, rho, g]).all()
            or time < 0 or water_depth <= 0 or ramp_time < 0
            or rho <= 0 or g <= 0 or spread[0] < 0):
        raise ValueError("source drag needs one finite irregular heading and body state")
    if not elements:
        raise ValueError("source drag needs Morison elements")
    k = finite_depth_wavenumber(omega, water_depth=water_depth, gravity=g)
    ramp = ((1 - np.cos(np.pi * time / ramp_time)) / 2
            if ramp_time and time < ramp_time else 1.)
    amplitude = np.sqrt(spread[0] * spectrum * width) * ramp
    roll, pitch, yaw = pose[3:]
    c4, s4 = np.cos(roll), np.sin(roll)
    c5, s5 = np.cos(pitch), np.sin(pitch)
    c6, s6 = np.cos(yaw), np.sin(yaw)
    rotation = np.array([
        [c5 * c6, c4 * s6 + s4 * s5 * c6, s4 * s6 - c4 * s5 * s6],
        [-c5 * s6, c4 * c6 - s4 * s5 * s6, s4 * c6 + c4 * s5 * s6],
        [s5, -s4 * c5, c4 * c5],
    ])
    heading_rad = np.deg2rad(heading[0])
    axis = np.array([np.cos(heading_rad), np.sin(heading_rad)])
    result = np.zeros(6)
    kh = k * water_depth
    for element in elements:
        point = np.asarray(element.point, dtype=float)
        cd = np.asarray(element.drag_coefficient, dtype=float)
        area = np.asarray(element.area, dtype=float)
        ca = np.asarray(element.added_mass_coefficient, dtype=float)
        if (any(value.shape != (3,) or not np.isfinite(value).all()
                for value in (point, cd, area, ca))
                or np.any(cd < 0) or np.any(area < 0)
                or np.any(ca != 0) or element.volume != 0):
            raise ValueError("source drag needs zero-volume, nonnegative elements")
        rotated = rotation @ point
        world = pose[:3] + rotated
        if world[2] > 0:
            continue
        kz = k * world[2]
        deep = kh > np.pi
        horizontal = np.where(deep, np.exp(kz),
                              np.cosh(kz + kh) / np.cosh(kh))
        vertical = np.where(deep, np.exp(kz),
                            np.sinh(kz + kh) / np.cosh(kh))
        argument = omega * time - k * (world[:2] @ axis) + phase[:, 0]
        fluid_horizontal = np.sum(
            amplitude * horizontal * np.cos(argument) * g * k / omega
        ) * axis
        fluid_vertical = np.sum(
            -amplitude * vertical * np.sin(argument) * g * k / omega
        )
        fluid_velocity = np.r_[fluid_horizontal, fluid_vertical]
        body_velocity = speed[:3] + np.cross(speed[3:], point)
        relative = fluid_velocity - body_velocity
        force = (0.5 * rho * (np.abs(cd) @ rotation)
                 * np.abs(area @ rotation) * relative * np.abs(relative))
        result[:3] += force
        result[3:] += np.cross(rotated, force)
    return result

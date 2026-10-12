"""Morison drag and added mass coupled to moving WEC bodies."""

from typing import Sequence

import numpy as np

from .core import (
    MorisonElement, finite_depth_wavenumber,
    _current_speed_at_depth, _validate_current,
)


def _axial_heave_element(element: MorisonElement):
    if not isinstance(element, MorisonElement):
        raise TypeError("elements must be MorisonElement values")
    point = np.asarray(element.point, dtype=float)
    cd = np.asarray(element.drag_coefficient, dtype=float)
    ca = np.asarray(element.added_mass_coefficient, dtype=float)
    area = np.asarray(element.area, dtype=float)
    if (any(value.shape != (3,) or not np.isfinite(value).all()
            for value in (point, cd, ca, area))
            or not np.allclose(point[:2], 0, rtol=0, atol=1e-12)
            or any(not np.allclose(value[:2], 0, rtol=0, atol=1e-12)
                   for value in (cd, ca, area))
            or np.any(cd < 0) or np.any(ca < 0) or np.any(area < 0)
            or not np.isfinite(element.volume) or element.volume <= 0
            or element.axis is not None):
        raise ValueError("heave Morison elements need axial, nonnegative geometry")
    return point[2], cd[2], ca[2], area[2], element.volume


def no_wave_heave_morison_terms(
    elements: Sequence[MorisonElement], *, center_z: float,
    heave: float, speed: float, rho: float,
) -> tuple[float, float]:
    """Return drag force and added mass for submerged axial heave elements.

    In still water, Cartesian option 1 reduces to ``-Cd*A*rho*v*|v|/2``
    plus ``-Ca*V*rho*a``. The latter belongs on the equation's mass side.
    Submergence follows WEC-Sim's mean-waterline switch at the element point.
    """
    if (not np.isfinite([center_z, heave, speed, rho]).all() or rho <= 0):
        raise ValueError("Morison heave state and density must be finite")
    drag_force = 0.0
    added_mass = 0.0
    for element in elements:
        point_z, cd, ca, area, volume = _axial_heave_element(element)
        if center_z + heave + point_z > 0:
            continue
        drag_force -= 0.5 * rho * cd * area * speed * abs(speed)
        added_mass += rho * volume * ca
    return drag_force, added_mass


def regular_wave_heave_morison_terms(
    elements: Sequence[MorisonElement], *, center_z: float,
    heave: float, speed: float, time: float, wave_height: float,
    wave_period: float, ramp_time: float, water_depth: float,
    rho: float, g: float = 9.81,
) -> tuple[float, float]:
    """Return axial regular-wave Morison force without body acceleration and mass.

    The element follows the moving heave coordinate. The fluid vertical
    velocity and acceleration are evaluated at its current submerged point.
    The returned added mass multiplies body acceleration in the dynamics.
    """
    if (not np.isfinite([center_z, heave, speed, time, wave_height,
                         wave_period, ramp_time, rho, g]).all()
            or time < 0 or wave_height < 0 or wave_period <= 0
            or ramp_time < 0 or rho <= 0 or g <= 0
            or not (np.isfinite(water_depth) or np.isposinf(water_depth))
            or water_depth <= 0):
        raise ValueError("regular-wave Morison heave settings are invalid")
    omega = 2 * np.pi / wave_period
    k = (omega * omega / g if np.isposinf(water_depth) else
         finite_depth_wavenumber(
             np.array([omega]), water_depth=water_depth, gravity=g,
         )[0])
    amplitude = wave_height / 2
    if ramp_time and time < ramp_time:
        amplitude *= (1 - np.cos(np.pi * time / ramp_time)) / 2
    force = 0.0
    added_mass = 0.0
    for element in elements:
        point_z, cd, ca, area, volume = _axial_heave_element(element)
        world_z = center_z + heave + point_z
        if world_z > 0:
            continue
        if k * water_depth > np.pi:
            vertical = np.exp(k * world_z)
        else:
            kh = k * water_depth
            vertical = np.sinh(k * world_z + kh) / np.cosh(kh)
        phase = omega * time
        fluid_speed = (-amplitude * vertical * np.sin(phase)
                       * g * k / omega)
        fluid_acceleration = (-amplitude * vertical * np.cos(phase)
                              * g * k)
        relative_speed = fluid_speed - speed
        force += (0.5 * rho * cd * area * relative_speed * abs(relative_speed)
                  + rho * volume * (1 + ca) * fluid_acceleration)
        added_mass += rho * volume * ca
    return force, added_mass


def regular_wave_axial_morison_terms(
    elements: Sequence[MorisonElement], *, position: Sequence[float],
    velocity: Sequence[float], time: float, wave_height: float,
    wave_period: float, ramp_time: float, water_depth: float,
    rho: float, g: float = 9.81,
    current_speed: float = 0.0, current_direction: float = 0.0,
    current_profile: str = "uniform", current_depth: float | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Physical axial-element wrench and added mass for surge/heave/pitch.

    The point and its axial direction rotate with body pitch. The returned
    wrench excludes the body's translational and angular acceleration; the
    six-by-six matrix places that acceleration dependence on the mass side.
    This uses a proper pitch rotation, independently of the pinned MATLAB
    source function's nonorthogonal general rotation diagnostic.
    """
    pose = np.asarray(position, dtype=float)
    speed = np.asarray(velocity, dtype=float)
    if (pose.shape != (6,) or speed.shape != (6,)
            or not np.isfinite(pose).all() or not np.isfinite(speed).all()
            or not np.allclose(pose[[1, 3, 5]], 0, rtol=0, atol=1e-12)
            or not np.allclose(speed[[1, 3, 5]], 0, rtol=0, atol=1e-12)):
        raise ValueError("axial moving Morison needs finite surge/heave/pitch state")
    if (not np.isfinite([time, wave_height, wave_period, ramp_time, rho, g]).all()
            or time < 0 or wave_height < 0 or wave_period <= 0
            or ramp_time < 0 or rho <= 0 or g <= 0
            or not (np.isfinite(water_depth) or np.isposinf(water_depth))
            or water_depth <= 0):
        raise ValueError("regular-wave Morison settings are invalid")
    _validate_current(current_speed, current_direction,
                      current_profile, current_depth)
    omega = 2 * np.pi / wave_period
    k = (omega * omega / g if np.isposinf(water_depth) else
         finite_depth_wavenumber(
             np.array([omega]), water_depth=water_depth, gravity=g,
         )[0])
    amplitude = wave_height / 2
    if ramp_time and time < ramp_time:
        amplitude *= (1 - np.cos(np.pi * time / ramp_time)) / 2
    pitch = pose[4]
    axis = np.array([np.sin(pitch), 0.0, np.cos(pitch)])
    angular_velocity = np.array([0.0, speed[4], 0.0])
    wrench = np.zeros(6)
    added_mass = np.zeros((6, 6))
    for element in elements:
        point_z, cd, ca, area, volume = _axial_heave_element(element)
        arm = point_z * axis
        world = pose[:3] + arm
        if world[2] > 0:
            continue
        kh = k * water_depth
        if kh > np.pi:
            horizontal = vertical = np.exp(k * world[2])
        else:
            kz = k * world[2]
            horizontal = np.cosh(kz + kh) / np.cosh(kh)
            vertical = np.sinh(kz + kh) / np.cosh(kh)
        phase = omega * time - k * world[0]
        wave_speed = amplitude * g * k / omega
        wave_acceleration = amplitude * g * k
        fluid_velocity = np.array([
            wave_speed * horizontal * np.cos(phase), 0.0,
            -wave_speed * vertical * np.sin(phase),
        ])
        fluid_velocity[:2] += _current_speed_at_depth(
            world[2], time, ramp_time, current_speed, current_profile,
            current_depth,
        ) * np.array([np.cos(np.deg2rad(current_direction)),
                      np.sin(np.deg2rad(current_direction))])
        fluid_acceleration = np.array([
            -wave_acceleration * horizontal * np.sin(phase), 0.0,
            -wave_acceleration * vertical * np.cos(phase),
        ])
        point_velocity = speed[:3] + np.cross(angular_velocity, arm)
        bias_acceleration = np.cross(
            angular_velocity, np.cross(angular_velocity, arm),
        )
        axial_speed = axis @ (fluid_velocity - point_velocity)
        force = (0.5 * rho * cd * area * axial_speed * abs(axial_speed) * axis
                 + rho * volume * fluid_acceleration
                 + rho * volume * ca * axis
                 * (axis @ (fluid_acceleration - bias_acceleration)))
        wrench[:3] += force
        wrench[3:] += np.cross(arm, force)
        added_mass[:3, :3] += rho * volume * ca * np.outer(axis, axis)
    return wrench, added_mass


def regular_wave_normal_tangential_morison_terms(
    elements: Sequence[MorisonElement], *, position: Sequence[float],
    velocity: Sequence[float], time: float, wave_height: float,
    wave_period: float, direction: float = 0.0, ramp_time: float,
    water_depth: float, rho: float, g: float = 9.81,
    current_speed: float = 0.0, current_direction: float = 0.0,
    current_profile: str = "uniform", current_depth: float | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Physical option-2 wrench and added mass on a six-DOF body state.

    The first two coefficient and area entries are normal and tangential;
    the third is unused, as in WEC-Sim option 2. The axis and body-local
    point rotate with a proper XYZ attitude matrix. The returned 6x6 mass includes
    point-offset translation/rotation coupling and is positive semidefinite
    for nonnegative added-mass coefficients. This is not the pinned source's
    nonorthogonal general rotation diagnostic.
    """
    pose = np.asarray(position, dtype=float)
    speed = np.asarray(velocity, dtype=float)
    if (pose.shape != (6,) or speed.shape != (6,)
            or not np.isfinite(pose).all() or not np.isfinite(speed).all()):
        raise ValueError("moving Morison needs finite six-DOF state")
    if (not np.isfinite([time, wave_height, wave_period, direction,
                         ramp_time, rho, g]).all()
            or time < 0 or wave_height < 0 or wave_period <= 0
            or ramp_time < 0 or rho <= 0 or g <= 0
            or not (np.isfinite(water_depth) or np.isposinf(water_depth))
            or water_depth <= 0 or not elements):
        raise ValueError("regular-wave Morison settings are invalid")
    _validate_current(current_speed, current_direction,
                      current_profile, current_depth)
    omega = 2 * np.pi / wave_period
    k = (omega * omega / g if np.isposinf(water_depth) else
         finite_depth_wavenumber(
             np.array([omega]), water_depth=water_depth, gravity=g,
         )[0])
    amplitude = wave_height / 2
    if ramp_time and time < ramp_time:
        amplitude *= (1 - np.cos(np.pi * time / ramp_time)) / 2
    roll, pitch, yaw = pose[3:]
    cr, sr = np.cos(roll), np.sin(roll)
    cp, sp = np.cos(pitch), np.sin(pitch)
    cy, sy = np.cos(yaw), np.sin(yaw)
    rotation = np.array([
        [cy * cp, cy * sp * sr - sy * cr, cy * sp * cr + sy * sr],
        [sy * cp, sy * sp * sr + cy * cr, sy * sp * cr - cy * sr],
        [-sp, cp * sr, cp * cr],
    ])
    angular_velocity = speed[3:]
    wave_axis = np.array([np.cos(np.deg2rad(direction)),
                          np.sin(np.deg2rad(direction))])
    wrench = np.zeros(6)
    added_mass = np.zeros((6, 6))
    for element in elements:
        if not isinstance(element, MorisonElement) or element.axis is None:
            raise ValueError("normal/tangential elements need a body-local axis")
        point = np.asarray(element.point, dtype=float)
        axis = np.asarray(element.axis, dtype=float)
        cd = np.asarray(element.drag_coefficient, dtype=float)
        ca = np.asarray(element.added_mass_coefficient, dtype=float)
        area = np.asarray(element.area, dtype=float)
        if (any(value.shape != (3,) or not np.isfinite(value).all()
                for value in (point, axis, cd, ca, area))
                or np.linalg.norm(axis) == 0
                or np.any(cd < 0) or np.any(ca < 0) or np.any(area < 0)
                or not np.isfinite(element.volume) or element.volume <= 0):
            raise ValueError("normal/tangential geometry and coefficients are invalid")
        arm = rotation @ point
        world = pose[:3] + arm
        if world[2] > 0:
            continue
        tangent = rotation @ axis
        tangent /= np.linalg.norm(tangent)
        tangential = np.outer(tangent, tangent)
        normal = np.eye(3) - tangential
        kh = k * water_depth
        if kh > np.pi:
            horizontal = vertical = np.exp(k * world[2])
        else:
            kz = k * world[2]
            horizontal = np.cosh(kz + kh) / np.cosh(kh)
            vertical = np.sinh(kz + kh) / np.cosh(kh)
        phase = omega * time - k * (world[:2] @ wave_axis)
        wave_speed = amplitude * g * k / omega
        wave_acceleration = amplitude * g * k
        fluid_velocity = np.r_[
            wave_speed * horizontal * np.cos(phase) * wave_axis,
            -wave_speed * vertical * np.sin(phase),
        ]
        fluid_velocity[:2] += _current_speed_at_depth(
            world[2], time, ramp_time, current_speed, current_profile,
            current_depth,
        ) * np.array([np.cos(np.deg2rad(current_direction)),
                      np.sin(np.deg2rad(current_direction))])
        fluid_acceleration = np.r_[
            -wave_acceleration * horizontal * np.sin(phase) * wave_axis,
            -wave_acceleration * vertical * np.cos(phase),
        ]
        point_velocity = speed[:3] + np.cross(angular_velocity, arm)
        bias_acceleration = np.cross(
            angular_velocity, np.cross(angular_velocity, arm),
        )
        relative = fluid_velocity - point_velocity
        relative_n = normal @ relative
        relative_t = tangential @ relative
        inertia = ca[0] * normal + ca[1] * tangential
        force = (0.5 * rho * (
            cd[0] * area[0] * relative_n * np.linalg.norm(relative_n)
            + cd[1] * area[1] * relative_t * np.linalg.norm(relative_t)
        ) + rho * element.volume * (
            fluid_acceleration
            + inertia @ (fluid_acceleration - bias_acceleration)
        ))
        wrench[:3] += force
        wrench[3:] += np.cross(arm, force)
        skew = np.array([[0, -arm[2], arm[1]],
                         [arm[2], 0, -arm[0]],
                         [-arm[1], arm[0], 0]])
        point_jacobian = np.column_stack((np.eye(3), -skew))
        added_mass += (rho * element.volume * point_jacobian.T
                       @ inertia @ point_jacobian)
    return wrench, added_mass

"""Axial cable tension law used by the published WEC-Sim Cable application."""

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class WecSimCableTension:
    """Reproduce ``calcCableTens`` in the cable's local z coordinate.

    ``relative_position`` and ``relative_velocity`` are the cable block's
    displacement and speed relative to its initial endpoint separation.
    The returned force uses the source block's signed actuation convention.
    The source does not clamp the damping term after extension, so a rapidly
    contracting stretched cable can report positive force.
    """

    stiffness: float
    damping: float
    length: float
    initial_length: float

    def __post_init__(self):
        parameters = np.asarray((
            self.stiffness, self.damping, self.length, self.initial_length,
        ), dtype=float)
        if (not np.isfinite(parameters).all()
                or np.any(parameters < 0)):
            raise ValueError("cable coefficients and lengths must be finite and nonnegative")

    def force_z(self, relative_position, relative_velocity):
        position = np.asarray(relative_position, dtype=float)
        velocity = np.asarray(relative_velocity, dtype=float)
        if not np.isfinite(position).all() or not np.isfinite(velocity).all():
            raise ValueError("cable position and velocity must be finite")
        stretched_length = np.abs(position + self.initial_length)
        force = np.where(
            stretched_length > self.length,
            -self.stiffness * (stretched_length - self.length)
            - self.damping * velocity,
            0.,
        )
        return float(force) if np.ndim(force) == 0 else force


@dataclass(frozen=True)
class PlanarCableAttachment:
    """Cable endpoints fixed to two body centers in the x/z pitch plane.

    Offsets are body-local ``(x, z)`` coordinates from each body's center of
    gravity. Poses and rates use ``(x, z, pitch)`` and ``(vx, vz, pitch_rate)``.
    The returned displacement is endpoint distance minus ``initial_length``.
    """

    base_offset: tuple[float, float]
    follower_offset: tuple[float, float]
    initial_length: float

    def __post_init__(self):
        for label in ("base_offset", "follower_offset"):
            value = np.asarray(getattr(self, label), dtype=float)
            if value.shape != (2,) or not np.isfinite(value).all():
                raise ValueError(f"{label} must be a finite x/z point")
        if not np.isfinite(self.initial_length) or self.initial_length <= 0:
            raise ValueError("initial cable length must be positive")

    def motion(self, base_pose, base_rate, follower_pose, follower_rate):
        arrays = [np.asarray(value, dtype=float) for value in (
            base_pose, base_rate, follower_pose, follower_rate,
        )]
        if (arrays[0].shape != (3,) and
                (arrays[0].ndim != 2 or arrays[0].shape[1] != 3)):
            raise ValueError("cable poses and rates must have three x/z/pitch columns")
        if any(value.shape != arrays[0].shape or not np.isfinite(value).all()
               for value in arrays):
            raise ValueError("cable poses and rates must be finite and aligned")

        def endpoint(pose, rate, offset):
            local_x, local_z = offset
            sine, cosine = np.sin(pose[..., 2]), np.cos(pose[..., 2])
            point = np.stack((
                pose[..., 0] + local_x * cosine + local_z * sine,
                pose[..., 1] - local_x * sine + local_z * cosine,
            ), axis=-1)
            speed = np.stack((
                rate[..., 0] + (-local_x * sine + local_z * cosine)
                * rate[..., 2],
                rate[..., 1] + (-local_x * cosine - local_z * sine)
                * rate[..., 2],
            ), axis=-1)
            return point, speed

        base_point, base_speed = endpoint(arrays[0], arrays[1],
                                          self.base_offset)
        follower_point, follower_speed = endpoint(arrays[2], arrays[3],
                                                  self.follower_offset)
        vector = follower_point - base_point
        length = np.linalg.norm(vector, axis=-1)
        if np.any(length <= 0):
            raise ValueError("cable endpoints must not coincide")
        speed = np.sum(vector * (follower_speed - base_speed), axis=-1) / length
        displacement = length - self.initial_length
        if np.ndim(displacement) == 0:
            return float(displacement), float(speed)
        return displacement, speed

    def world_wrenches(self, base_pose, follower_pose, force_z):
        """Map signed axial cable force to body-center ``(Fx, Fz, My)``.

        Positive ``force_z`` acts on the follower from base toward follower.
        The base receives the equal and opposite force. Moments use each
        body's local attachment arm and WEC-Sim's positive pitch convention.
        """
        base, follower = (np.asarray(pose, dtype=float)
                          for pose in (base_pose, follower_pose))
        force = np.asarray(force_z, dtype=float)
        if (base.shape != follower.shape
                or (base.shape != (3,) and
                    (base.ndim != 2 or base.shape[1] != 3))
                or force.shape != base.shape[:-1]
                or not all(np.isfinite(value).all()
                           for value in (base, follower, force))):
            raise ValueError("cable poses and force must be finite and aligned")

        def endpoint(pose, offset):
            local_x, local_z = offset
            sine, cosine = np.sin(pose[..., 2]), np.cos(pose[..., 2])
            arm = np.stack((local_x * cosine + local_z * sine,
                            -local_x * sine + local_z * cosine), axis=-1)
            return pose[..., :2] + arm, arm

        base_point, base_arm = endpoint(base, self.base_offset)
        follower_point, follower_arm = endpoint(follower, self.follower_offset)
        direction = follower_point - base_point
        distance = np.linalg.norm(direction, axis=-1)
        if np.any(distance <= 0):
            raise ValueError("cable endpoints must not coincide")
        follower_force = force[..., None] * direction / distance[..., None]
        base_force = -follower_force

        def wrench(applied, arm):
            moment_y = arm[..., 1] * applied[..., 0] - arm[..., 0] * applied[..., 1]
            return np.concatenate((applied, moment_y[..., None]), axis=-1)

        return wrench(base_force, base_arm), wrench(follower_force, follower_arm)

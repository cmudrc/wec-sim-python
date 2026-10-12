"""Smooth rotational joint friction used by PTO and linkage models."""

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class StribeckFriction:
    """Breakaway, Coulomb, and viscous torque opposing angular motion."""

    breakaway_torque: float
    breakaway_speed: float
    coulomb_torque: float
    viscous_damping: float

    def __post_init__(self):
        values = (self.breakaway_torque, self.breakaway_speed,
                  self.coulomb_torque, self.viscous_damping)
        if (not np.isfinite(values).all() or self.breakaway_speed <= 0
                or self.coulomb_torque < 0 or self.viscous_damping < 0
                or self.breakaway_torque < self.coulomb_torque):
            raise ValueError("rotational friction settings are invalid")

    def torque(self, angular_speed):
        """Return passive torque, in N m, at speed in rad/s."""
        speed = -np.asarray(angular_speed, dtype=float)
        static_scale = np.sqrt(2 * np.e) * (
            self.breakaway_torque - self.coulomb_torque
        )
        static_threshold = np.sqrt(2) * self.breakaway_speed
        coulomb_threshold = self.breakaway_speed / 10
        ratio = speed / static_threshold
        return (self.viscous_damping * speed
                + static_scale * ratio * np.exp(-ratio**2)
                + self.coulomb_torque * np.tanh(speed / coulomb_threshold))

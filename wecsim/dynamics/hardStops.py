"""Unilateral spring-damper limits for a translational PTO stroke."""

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class LinearHardStops:
    """Force that pushes a PTO stroke back between its lower and upper bounds.

    Each active stop has a linear spring beyond its bound. Within the narrow
    transition width, a cubic smoothstep raises the force continuously. The
    contact damper cannot pull the slider farther into the stop on rebound.
    """

    lower_bound: float
    upper_bound: float
    lower_stiffness: float
    upper_stiffness: float
    lower_damping: float = 0.0
    upper_damping: float = 0.0
    lower_transition_width: float = 1e-4
    upper_transition_width: float = 1e-4

    def __post_init__(self):
        values = np.asarray(tuple(self.__dict__.values()), dtype=float)
        if not np.isfinite(values).all():
            raise ValueError("hard-stop settings must be finite")
        if self.lower_bound >= self.upper_bound:
            raise ValueError("hard-stop lower bound must precede upper bound")
        if self.lower_stiffness <= 0 or self.upper_stiffness <= 0:
            raise ValueError("hard-stop stiffness must be positive")
        if self.lower_damping < 0 or self.upper_damping < 0:
            raise ValueError("hard-stop damping must be nonnegative")
        if self.lower_transition_width <= 0 or self.upper_transition_width <= 0:
            raise ValueError("hard-stop transition width must be positive")

    @staticmethod
    def _smoothing(penetration, width):
        fraction = np.clip(penetration / width, 0.0, 1.0)
        return fraction * fraction * (3 - 2 * fraction)

    def force(self, stroke, speed):
        """Return signed stop force in the positive-stroke direction (N)."""
        stroke, speed = np.broadcast_arrays(
            np.asarray(stroke, dtype=float), np.asarray(speed, dtype=float),
        )
        upper = np.maximum(stroke - self.upper_bound, 0.0)
        lower = np.maximum(self.lower_bound - stroke, 0.0)
        upper_force = self._smoothing(upper, self.upper_transition_width) * (
            self.upper_stiffness * upper + self.upper_damping * speed
        )
        lower_force = self._smoothing(lower, self.lower_transition_width) * (
            self.lower_stiffness * lower - self.lower_damping * speed
        )
        return np.maximum(lower_force, 0.0) - np.maximum(upper_force, 0.0)

"""Rotary-to-linear linkages in published OSWEC applications."""

from dataclasses import dataclass

import numpy as np


def _initial_horizontal(crank: float, offset: float, rod_length: float) -> float:
    parameters = np.asarray((crank, offset, rod_length), dtype=float)
    if (not np.isfinite(parameters).all() or crank <= 0 or rod_length <= 0
            or rod_length <= abs(crank - offset)):
        raise ValueError("crank, offset, and initial rod geometry are invalid")
    return float(np.sqrt(rod_length**2 - (crank - offset)**2))


@dataclass(frozen=True)
class FixedRodCrank:
    """Published slider-crank position measured from its zero-angle pose."""

    crank: float
    offset: float
    rod_length: float

    def __post_init__(self):
        _initial_horizontal(self.crank, self.offset, self.rod_length)

    def stroke_and_jacobian(self, angle):
        theta = np.asarray(angle, dtype=float)
        horizontal = self.crank * np.cos(theta) - self.offset
        under_root = self.rod_length**2 - horizontal**2
        if not np.isfinite(theta).all() or np.any(under_root <= 0):
            raise ValueError("crank angle is outside the slider geometry")
        rod_horizontal = np.sqrt(under_root)
        stroke = (self.crank * np.sin(theta) + rod_horizontal
                  - _initial_horizontal(self.crank, self.offset,
                                        self.rod_length))
        jacobian = (self.crank * np.cos(theta)
                    + self.crank * horizontal * np.sin(theta)
                    / rod_horizontal)
        return stroke, jacobian

    def torque(self, piston_force, angle):
        """Published fixed-crank PTO torque from signed cylinder force."""
        return -np.asarray(piston_force, dtype=float) * self.stroke_and_jacobian(angle)[1]


@dataclass(frozen=True)
class AdjustableRodCrank:
    """Published adjustable-rod length change from its zero-angle length."""

    crank: float
    offset: float
    initial_rod_length: float

    def __post_init__(self):
        _initial_horizontal(self.crank, self.offset,
                            self.initial_rod_length)

    def stroke_and_jacobian(self, angle):
        theta = np.asarray(angle, dtype=float)
        if not np.isfinite(theta).all():
            raise ValueError("crank angle must be finite")
        initial_horizontal = _initial_horizontal(
            self.crank, self.offset, self.initial_rod_length,
        )
        dx = self.crank * np.cos(theta) - self.offset
        dy = initial_horizontal - self.crank * np.sin(theta)
        length = np.hypot(dx, dy)
        if np.any(length <= 0):
            raise ValueError("adjustable rod collapsed to zero length")
        stroke = length - self.initial_rod_length
        jacobian = -self.crank * (dx * np.sin(theta)
                                  + dy * np.cos(theta)) / length
        return stroke, jacobian

    def torque(self, piston_force, angle):
        """Published adjustable-rod PTO torque from signed cylinder force."""
        return np.asarray(piston_force, dtype=float) * self.stroke_and_jacobian(angle)[1]


@dataclass(frozen=True)
class PitchRodLinkage:
    """Straight rod from a body-local point to a fixed world anchor, in x/z.

    The body center, hinge, and anchor use world coordinates at zero pitch;
    ``body_point`` is relative to the body center. Positive pitch rotates
    positive z toward positive x, as in the OSWEC hinge cases.
    """

    anchor: tuple[float, float]
    hinge: tuple[float, float]
    body_center: tuple[float, float]
    body_point: tuple[float, float]

    def __post_init__(self):
        points = np.asarray((self.anchor, self.hinge, self.body_center,
                             self.body_point), dtype=float)
        if points.shape != (4, 2) or not np.isfinite(points).all():
            raise ValueError("linkage points must be finite x/z pairs")
        if np.linalg.norm(points[0] - points[2] - points[3]) <= 0:
            raise ValueError("rod length at zero pitch must be positive")

    def stroke_and_jacobian(self, angle):
        theta = np.asarray(angle, dtype=float)
        if not np.isfinite(theta).all():
            raise ValueError("pitch angle must be finite")
        anchor = np.asarray(self.anchor)
        hinge = np.asarray(self.hinge)
        offset = (np.asarray(self.body_center)
                  + np.asarray(self.body_point) - hinge)
        cosine, sine = np.cos(theta), np.sin(theta)
        dx = anchor[0] - hinge[0] - cosine * offset[0] - sine * offset[1]
        dz = anchor[1] - hinge[1] + sine * offset[0] - cosine * offset[1]
        length = np.hypot(dx, dz)
        if np.any(length <= 0):
            raise ValueError("rod collapsed to zero length")
        dx_dtheta = sine * offset[0] - cosine * offset[1]
        dz_dtheta = cosine * offset[0] + sine * offset[1]
        reference = np.linalg.norm(
            anchor - np.asarray(self.body_center) - np.asarray(self.body_point)
        )
        return length - reference, (dx * dx_dtheta + dz * dz_dtheta) / length

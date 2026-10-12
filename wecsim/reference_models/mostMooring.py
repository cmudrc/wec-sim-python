"""Static catenary restoring loads for the published MOST VolturnUS case."""

from dataclasses import dataclass

import numpy as np
from scipy.optimize import least_squares, root


@dataclass(frozen=True)
class MostStaticMooring:
    """Three equally spaced catenary lines with body-local fairleads.

    Defaults reproduce the pinned MOST input. ``pose`` is surge, sway,
    heave, roll, pitch, yaw relative to the mooring reference, in m and rad.
    """

    fairlead: tuple[float, float, float] = (-58.0, 0.0, -14.0)
    anchor: tuple[float, float, float] = (-837.0, 0.0, -200.0)
    line_length: float = 850.0
    diameter: float = .333
    mass_per_length: float = 685.0
    axial_stiffness: float = 3.27e9
    seabed_friction: float = 1.0
    water_density: float = 1025.0
    gravity: float = 9.80665

    def _line_geometry(self, tension, height, distance):
        horizontal, vertical = tension
        weight = (self.mass_per_length
                  - np.pi * self.diameter**2 / 4 * self.water_density) * self.gravity
        bottom_length = self.line_length - vertical / weight
        if bottom_length > 0:
            g = bottom_length - horizontal / (self.seabed_friction * weight)
            positive_g = max(g, 0)
            x = (bottom_length
                 + horizontal / weight
                 * np.arcsinh(weight * (self.line_length-bottom_length) / horizontal)
                 + horizontal * self.line_length / self.axial_stiffness
                 + self.seabed_friction * weight / (2*self.axial_stiffness)
                 * (g*positive_g-bottom_length**2))
            z = (horizontal / weight
                 * (np.sqrt(1+(weight*(self.line_length-bottom_length)/horizontal)**2)-1)
                 + weight*(self.line_length-bottom_length)**2
                 / (2*self.axial_stiffness))
        else:
            anchor_vertical = vertical-weight*self.line_length
            x = (horizontal / weight
                 * (np.arcsinh(vertical/horizontal)
                    - np.arcsinh(anchor_vertical/horizontal))
                 + horizontal*self.line_length/self.axial_stiffness)
            z = (horizontal / weight
                 * (np.sqrt(1+(vertical/horizontal)**2)
                    - np.sqrt(1+(anchor_vertical/horizontal)**2))
                 + (anchor_vertical*self.line_length
                    + weight*self.line_length**2/2)/self.axial_stiffness)
        return np.array([x-distance, z-height])

    @staticmethod
    def _rotation(roll, pitch, yaw):
        cx, sx = np.cos(roll), np.sin(roll)
        cy, sy = np.cos(pitch), np.sin(pitch)
        cz, sz = np.cos(yaw), np.sin(yaw)
        return np.array([
            [cy*cz, sx*sy*cz-cx*sz, cx*sy*cz+sx*sz],
            [cy*sz, sx*sy*sz+cx*cz, cx*sy*sz-sx*cz],
            [-sy, sx*cy, cx*cy],
        ])

    def force(self, pose) -> tuple[np.ndarray, np.ndarray]:
        """Return six-component platform load and each line's H/V tension."""
        pose = np.asarray(pose, dtype=float)
        if pose.shape != (6,) or not np.isfinite(pose).all():
            raise ValueError("pose must contain six finite coordinates")
        rotation = self._rotation(*pose[3:])
        total = np.zeros(6)
        tensions = np.zeros((3, 2))
        for line, beta in enumerate(2*np.pi*np.arange(3)/3):
            cb, sb = np.cos(beta), np.sin(beta)
            azimuth = np.array([[cb, -sb, 0], [sb, cb, 0], [0, 0, 1]])
            fairlead = rotation @ (azimuth @ np.asarray(self.fairlead))
            anchor = azimuth @ np.asarray(self.anchor)
            separation = anchor - fairlead - pose[:3]
            height = abs(separation[2])
            distance = np.linalg.norm(separation[:2])
            solution = root(
                self._line_geometry, [1e6, 2e6],
                args=(height, distance), method="hybr",
                options={"xtol": 1e-12},
            )
            # Fall back when the fast solve misses a positive, accurate root.
            if (not np.isfinite(solution.x).all()
                    or not (solution.x > 1e-9).all()
                    or not np.isfinite(solution.fun).all()
                    or np.linalg.norm(solution.fun) > 1e-6):
                solution = least_squares(
                    self._line_geometry, [1e6, 2e6],
                    args=(height, distance), bounds=(1e-9, np.inf),
                    xtol=1e-12, ftol=1e-12, gtol=1e-12, max_nfev=300,
                )
            if np.linalg.norm(solution.fun) > 1e-6:
                raise RuntimeError("static mooring geometry did not converge")
            horizontal, vertical = solution.x
            direction = separation[:2] / distance
            force = np.array([horizontal*direction[0],
                              horizontal*direction[1], -vertical])
            total += np.r_[force, np.cross(fairlead, force)]
            tensions[line] = solution.x
        return total, tensions

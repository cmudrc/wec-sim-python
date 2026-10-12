"""Read the TurbSim full-field wind input used by the MOST application."""

from dataclasses import dataclass
from math import ceil, floor
from pathlib import Path
import struct

import numpy as np


_HEADER = struct.Struct("<h4i12fi")


@dataclass(frozen=True)
class TurbSimWind:
    """Decoded wind; velocity axes are time, component (U/V/W), y, z."""

    velocity: np.ndarray
    tower_velocity: np.ndarray
    y: np.ndarray
    z: np.ndarray
    dt: float
    hub_height: float
    mean_wind_speed: float
    description: str


@dataclass(frozen=True)
class MostConstantWind:
    """Uniform wind over MOST's published constant-wind domain."""

    speed: float
    end_time: float
    direction: tuple[float, float, float] = (1.0, 0.0, 0.0)
    domain_sizes: tuple[float, float, float] = (1e4, 1e4, 1e4)

    def __post_init__(self):
        direction = np.asarray(self.direction, dtype=float)
        domain = np.asarray(self.domain_sizes, dtype=float)
        if (not np.isfinite([self.speed, self.end_time]).all()
                or self.speed < 0 or self.end_time <= 0
                or direction.shape != (3,) or not np.isfinite(direction).all()
                or np.linalg.norm(direction) == 0
                or domain.shape != (3,) or not np.isfinite(domain).all()
                or np.any(domain <= 0)):
            raise ValueError("MOST constant wind needs finite speed, time, direction, and domain")

    def sampler_at(self, time: float):
        """Return the spatially uniform vector at a valid simulation time."""
        if not np.isfinite(time) or not 0 <= time <= self.end_time:
            raise ValueError("MOST wind time is outside the available record")
        direction = np.asarray(self.direction, dtype=float)
        vector = self.speed * direction / np.linalg.norm(direction)
        domain = np.asarray(self.domain_sizes, dtype=float)

        def sample(position):
            point = np.asarray(position, dtype=float)
            if point.shape != (3,) or not np.isfinite(point).all():
                raise ValueError("MOST wind position must be a finite 3-vector")
            if (abs(point[0]) > domain[0] / 2
                    or abs(point[1]) > domain[1] / 2
                    or not 0 <= point[2] <= domain[2]):
                return np.full(3, np.nan)
            return vector.copy()

        return sample


class MostWindField:
    """MOST's frozen-turbulence X planes, evaluated without a large 5-D copy.

    ``at_index(i)`` has axes component, X, Y, Z, matching one time step of
    MOST's ``Wind.SpatialDiscrUVW`` from ``RunTurbsim.m``.
    """

    def __init__(self, wind: TurbSimWind, speed: float = 8.0):
        if not np.isfinite(speed) or speed <= 0:
            raise ValueError("wind advection speed must be positive and finite")
        y = np.asarray(wind.y, dtype=float)
        z = np.asarray(wind.z, dtype=float)
        if (y.ndim != 1 or z.ndim != 1 or not len(y) or not len(z)
                or not np.isfinite(y).all() or not np.isfinite(z).all()
                or np.any(np.diff(y) <= 0) or np.any(np.diff(z) <= 0)):
            raise ValueError("MOST wind coordinates must increase")
        self.wind = wind
        self.speed = float(speed)
        self.x = np.arange(-30.0, 61.0, 10.0)
        self.discarded = ceil(wind.y[-1] / (self.speed * wind.dt))
        self.n_time = len(wind.velocity) - 2 * self.discarded + 1
        if self.n_time < 1:
            raise ValueError("TurbSim record is too short for MOST advection")
        first = floor(self.discarded + 1 - self.x[-1] / (self.speed * wind.dt) + .5)
        last = floor(self.discarded + self.n_time - self.x[0] / (self.speed * wind.dt) + .5)
        if first < 1 or last > len(wind.velocity):
            raise ValueError("MOST advection samples outside the TurbSim record")

    @property
    def time(self) -> np.ndarray:
        return np.arange(self.n_time) * self.wind.dt

    def at_index(self, index: int) -> np.ndarray:
        """Return all ten X planes at a zero-based output time index."""
        if not 0 <= index < self.n_time:
            raise IndexError("MOST wind time index is out of range")
        source = np.floor(
            self.discarded + 1 + index - self.x / (self.speed * self.wind.dt) + .5,
        ).astype(int) - 1
        return self.wind.velocity[source].transpose(1, 0, 2, 3)

    def sampler(self, index: int):
        """Return MOST's trilinear world-position wind sampler at one time.

        Outside the published X/Y/Z grid the sampler returns NaN, as the
        MATLAB BEM block's ``interp3`` call does.
        """
        values = self.at_index(index).transpose(1, 2, 3, 0)
        return self._sampler_from_values(values)

    def _sampler_from_values(self, values):
        """Interpolate one blade-element position without grid-object overhead."""
        grids = (self.x, self.wind.y, self.wind.z)

        def sample(position):
            point = np.asarray(position, dtype=float)
            if point.shape != (3,):
                raise ValueError("MOST wind position must be a three-vector")
            if (not np.isfinite(point).all()
                    or any(point[axis] < grid[0] or point[axis] > grid[-1]
                           for axis, grid in enumerate(grids))):
                return np.full(3, np.nan)
            indices = []
            fractions = []
            for axis, grid in enumerate(grids):
                if len(grid) == 1:
                    indices.append((0, 0))
                    fractions.append(0.0)
                    continue
                lower = min(int(np.searchsorted(grid, point[axis], side="right")) - 1,
                            len(grid) - 2)
                indices.append((lower, lower + 1))
                fractions.append((point[axis] - grid[lower])
                                 / (grid[lower + 1] - grid[lower]))
            (x0, x1), (y0, y1), (z0, z1) = indices
            x, y, z = fractions
            return (
                (1-x)*(1-y)*(1-z)*values[x0, y0, z0]
                + (1-x)*(1-y)*z*values[x0, y0, z1]
                + (1-x)*y*(1-z)*values[x0, y1, z0]
                + (1-x)*y*z*values[x0, y1, z1]
                + x*(1-y)*(1-z)*values[x1, y0, z0]
                + x*(1-y)*z*values[x1, y0, z1]
                + x*y*(1-z)*values[x1, y1, z0]
                + x*y*z*values[x1, y1, z1]
            )

        return sample

    def sampler_at(self, time: float):
        """Sample the spatial field at a simulation time in seconds.

        MOST's time-series wind input is linearly interpolated between the
        TurbSim advection frames before the BEM spatial interpolation. Its
        Simulink input holds the final frame after the last recorded time.
        """
        if not np.isfinite(time) or time < 0:
            raise ValueError("MOST wind time is outside the available record")
        if time >= (self.n_time - 1) * self.wind.dt:
            return self.sampler(self.n_time - 1)
        fractional = time / self.wind.dt
        before = min(int(np.floor(fractional)), self.n_time - 1)
        after = min(before + 1, self.n_time - 1)
        weight = fractional - before
        if after == before or weight == 0:
            return self.sampler(before)
        values = ((1 - weight) * self.at_index(before)
                  + weight * self.at_index(after)).transpose(1, 2, 3, 0)
        return self._sampler_from_values(values)


def read_turbsim_bts(path: str | Path) -> TurbSimWind:
    """Read a little-endian TurbSim int16 ``.bts`` file.

    The returned grid ordering and scaling follow MOST's ``readfile_BTS.m``.
    Tower points are decoded even though the MOST reader does not return them.
    """
    data = Path(path).read_bytes()
    if len(data) < _HEADER.size:
        raise ValueError("TurbSim file has an incomplete header")
    (format_id, nz, ny, ntwr, nt, dz, dy, dt, mean_speed, hub_height, z1,
     slope_u, offset_u, slope_v, offset_v, slope_w, offset_w,
     description_length) = _HEADER.unpack_from(data)
    if format_id not in (7, 8):
        raise ValueError(f"unsupported TurbSim format identifier {format_id}")
    if min(nz, ny, nt) <= 0 or ntwr < 0 or not 0 <= description_length <= 4096:
        raise ValueError("invalid TurbSim grid dimensions or description length")
    slopes = np.array([slope_u, slope_v, slope_w], dtype=float)
    offsets = np.array([offset_u, offset_v, offset_w], dtype=float)
    if not np.isfinite(slopes).all() or np.any(slopes == 0):
        raise ValueError("invalid TurbSim velocity slopes")
    payload_start = _HEADER.size + description_length
    points_per_step = 3 * (ny * nz + ntwr)
    expected_size = payload_start + 2 * nt * points_per_step
    if len(data) != expected_size:
        raise ValueError(f"TurbSim file size is {len(data)}, expected {expected_size}")
    description = data[_HEADER.size:payload_start].decode("ascii", errors="replace")
    raw = np.frombuffer(data, dtype="<i2", count=nt * points_per_step,
                        offset=payload_start).reshape(nt, points_per_step)
    grid = raw[:, :3 * ny * nz].reshape(nt, nz, ny, 3).transpose(0, 3, 2, 1)
    velocity = (grid - offsets[None, :, None, None]) / slopes[None, :, None, None]
    tower = raw[:, 3 * ny * nz:].reshape(nt, ntwr, 3).transpose(0, 2, 1)
    tower_velocity = (tower - offsets[None, :, None]) / slopes[None, :, None]
    y = np.arange(ny) * dy - dy * (ny - 1) / 2
    z = np.arange(nz) * dz + z1
    return TurbSimWind(velocity, tower_velocity, y, z, dt, hub_height,
                       mean_speed, description)

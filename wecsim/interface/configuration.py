"""WEC body, wave, joint, and PTO specifications."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Sequence

import numpy as np

from wecsim.components.directLinearGenerator import DirectLinearGenerator
from wecsim.components.floatingOwc import FloatingOwcChamber, FloatingOwcTurbine
from wecsim.control.controls import DeclutchingControl, LatchingControl
from wecsim.dynamics.hardStops import LinearHardStops
from wecsim.dynamics.moorDyn import MoorDyn
from wecsim.hydrodynamics.morison import finite_depth_wavenumber


@dataclass(frozen=True)
class WorldPoint:
    """A fixed point in world coordinates, also usable as a rotation pivot."""

    x: float
    y: float
    z: float

    def coordinates(self) -> list[float]:
        return [self.x, self.y, self.z]


@dataclass(frozen=True)
class BodyPoint:
    """A point in a body's reference frame, relative to its center of gravity."""

    body: Body
    x: float
    y: float
    z: float

    def coordinates(self) -> list[float]:
        return [self.x, self.y, self.z]


@dataclass(frozen=True)
class Motion:
    body: Body
    dof: str
    scale: float = 1.0
    pivot: WorldPoint | BodyPoint | None = None


@dataclass(frozen=True)
class HydroState:
    """One BEM dataset and rigid properties for a variable-draft body."""

    hydro_file: str | Path
    mass: str | float = "equilibrium"
    inertia: tuple[float, float, float] = (0, 0, 0)


@dataclass(frozen=True)
class VariableHydro:
    """Ordered states; each switch time activates the next state."""

    states: tuple[HydroState, ...]
    switch_times: tuple[float, ...]


@dataclass(frozen=True)
class Body:
    name: str
    hydro_file: str | Path | None
    mass: str | float = "equilibrium"
    inertia: tuple[float, float, float] = (0, 0, 0)
    hydro_body: int | None = None
    mean_drift: str = "none"
    passive_yaw: bool = False
    geometry_file: str | Path | None = None
    nonlinear_hydro: str | None = None
    drag_coefficient: float = 0.0
    drag_area: float = 0.0
    variable_hydro: VariableHydro | None = None
    passive_yaw_threshold: float = 0.0
    yaw_heading_bank: tuple[float, ...] | None = None
    fixed: bool = False
    center_gravity: tuple[float, float, float] | None = None
    volume: float = 0.0

    def at(self, x: float, y: float, z: float) -> BodyPoint:
        """Locate an attachment or Morison point relative to this body's CG."""
        return BodyPoint(self, x, y, z)

    def move(self, dof: str, *, scale: float = 1.0,
             pivot: WorldPoint | BodyPoint | None = None) -> Motion:
        """Describe this body's contribution to one device coordinate."""
        return Motion(self, dof, scale, pivot)


@dataclass(frozen=True)
class _FloatingJoint:
    float_body: Body
    spar_body: Body
    location: WorldPoint
    pto_name: str
    damping: float
    stiffness: float
    equilibrium_position: float
    mooring_surge_stiffness: float
    hard_stops: LinearHardStops | None
    radiation_method: str | None
    added_mass_scheme: str
    moordyn: MoorDyn | None
    moordyn_point: BodyPoint | None


@dataclass(frozen=True)
class _FixedHinge:
    flap: Body
    base: Body | None
    location: WorldPoint
    pto_location: WorldPoint
    pto_name: str
    damping: float
    stiffness: float
    equilibrium_angle: float
    added_mass_scheme: str


@dataclass(frozen=True)
class _FloatingOwc:
    floater: Body
    column: Body
    moordyn: MoorDyn
    moordyn_point: BodyPoint | None
    column_height: float
    column_diameter: float
    pto_name: str
    pto_stiffness: float
    pto_damping: float
    chamber: FloatingOwcChamber | None
    turbine: FloatingOwcTurbine | None
    initial_rotor_speed: float


@dataclass(frozen=True)
class Coordinate:
    name: str
    motions: tuple[Motion, ...]


@dataclass(frozen=True)
class LinearPTO:
    name: str
    from_point: BodyPoint | WorldPoint
    to_point: BodyPoint | WorldPoint
    axis: tuple[float, float, float] | None = None
    damping: float = 0.0
    stiffness: float = 0.0
    equilibrium_position: float | None = None
    pretension: float | None = None
    control: DeclutchingControl | LatchingControl | None = None
    direct_drive: SimpleDirectDrive | None = None
    linear_generator: DirectLinearGenerator | None = None
    velocity_control_gain: float = 0.0


@dataclass(frozen=True)
class SimpleDirectDrive:
    """Reactive PI control and simple generator/drivetrain parameters."""

    kp: float
    ki: float
    torque_constant: float
    gear_ratio: float
    drivetrain_inertia: float
    drivetrain_friction: float
    winding_resistance: float
    winding_inductance: float


@dataclass(frozen=True)
class RotationalPTO:
    """A torsional spring and damper acting on a rotation coordinate."""

    name: str
    coordinate: Coordinate
    damping: float = 0.0
    stiffness: float = 0.0
    equilibrium_angle: float = 0.0


@dataclass(frozen=True)
class RegularWave:
    height: float
    period: float
    direction: float = 0.0
    current: Current | None = None

    def as_case(self) -> dict:
        wave = {"type": "regular", "height": self.height,
                "period": self.period, "direction": self.direction}
        if self.current is not None:
            if not isinstance(self.current, Current):
                raise TypeError("RegularWave.current must be Current")
            wave["current"] = self.current.as_case()
        return wave

    def elevation_at(
        self, time: Sequence[float], locations: Sequence[Sequence[float]], *,
        water_depth: float = np.inf, ramp_time: float = 0.0,
        g: float = 9.81,
    ) -> np.ndarray:
        """Sample the undisturbed, ramped surface at world XY locations.

        Rows follow ``time`` and columns follow ``locations``. This is the
        wave-marker signal; it does not depend on the WEC body motion.
        """
        times = np.asarray(time, dtype=float)
        points = np.asarray(locations, dtype=float)
        if (times.ndim != 1 or not np.isfinite(times).all()
                or np.any(times < 0) or points.ndim != 2
                or points.shape[1] != 2 or not np.isfinite(points).all()):
            raise ValueError("time and marker locations need finite 1D and Nx2 arrays")
        if (not np.isfinite([self.height, self.period, self.direction,
                             ramp_time, g]).all()
                or self.height < 0 or self.period <= 0 or ramp_time < 0
                or g <= 0 or np.isnan(water_depth) or water_depth <= 0):
            raise ValueError("regular wave and ramp settings must be valid")
        omega = 2 * np.pi / self.period
        k = (omega**2 / g if np.isinf(water_depth) else
             finite_depth_wavenumber(np.array([omega]),
                                     water_depth=water_depth, gravity=g)[0])
        heading = np.deg2rad(self.direction)
        distance = points[:, 0] * np.cos(heading) + points[:, 1] * np.sin(heading)
        ramp = np.ones_like(times)
        if ramp_time:
            starting = times < ramp_time
            ramp[starting] = 0.5 * (1 - np.cos(np.pi * times[starting] / ramp_time))
        return (self.height / 2 * ramp[:, None]
                * np.cos(omega * times[:, None] - k * distance[None, :]))


@dataclass(frozen=True)
class RegularCICWave:
    """Regular incident waves with convolution-integral radiation."""

    height: float
    period: float
    direction: float = 0.0

    def as_case(self) -> dict:
        return {"type": "regularCIC", "height": self.height,
                "period": self.period, "direction": self.direction}

    def elevation_at(self, time: Sequence[float],
                     locations: Sequence[Sequence[float]], *,
                     water_depth: float = np.inf, ramp_time: float = 0.0,
                     g: float = 9.81) -> np.ndarray:
        return RegularWave(self.height, self.period, self.direction).elevation_at(
            time, locations, water_depth=water_depth, ramp_time=ramp_time, g=g,
        )


@dataclass(frozen=True)
class Current:
    """Horizontal current for supported fixed or moving Morison bodies.

    ``profile`` is ``uniform``, ``power`` (the 1/7 law), or ``linear``.
    The latter two profiles need the current depth in metres.
    """

    speed: float
    direction: float = 0.0
    profile: str = "uniform"
    depth: float | None = None

    def as_case(self) -> dict:
        try:
            speed = float(self.speed)
            direction = float(self.direction)
            depth = None if self.depth is None else float(self.depth)
        except (TypeError, ValueError) as exc:
            raise ValueError("current needs numeric speed, direction, and depth") from exc
        if (not np.isfinite([speed, direction]).all()
                or speed < 0 or not -360 <= direction <= 360
                or self.profile not in ("uniform", "power", "linear")
                or (depth is not None and (not np.isfinite(depth) or depth <= 0))
                or (self.profile != "uniform" and depth is None)):
            raise ValueError("current needs valid speed, direction, profile, and depth")
        result = {"speed": speed, "direction": direction,
                  "profile": self.profile}
        if depth is not None:
            result["depth"] = depth
        return result


@dataclass(frozen=True)
class PMWave:
    """Pierson–Moskowitz sea for supported hydrodynamic or Morison bodies.

    ``height`` is significant wave height in metres, ``period`` is peak
    period in seconds, and ``direction`` is the incident heading in degrees.
    A saved phase CSV replays a MATLAB realization. By default, ``seed``
    selects a reproducible NumPy realization. Set ``phase_generator="matlab"``
    to use the pinned WEC-Sim Threefry substream with that seed.
    ``discretization="traditional"`` selects WEC-Sim's uniform frequency
    grid (1000 bins by default); the equal-energy default uses 500 bins.
    ``frequency_range`` narrows the BEM frequency interval in rad/s. A
    ``water_depth`` override is currently used only for fixed Morison bodies.
    """

    height: float
    period: float
    direction: float = 0.0
    seed: int | None = None
    phase_file: str | Path | None = None
    frequency_count: int | None = None
    directions: tuple[float, ...] | None = None
    spreading: tuple[float, ...] | None = None
    frequency_range: tuple[float, float] | None = None
    water_depth: float | None = None
    current: Current | None = None
    phase_generator: str = field(default="numpy", kw_only=True)
    discretization: str = field(default="equal_energy", kw_only=True)

    def as_case(self) -> dict:
        if self.seed is not None and self.phase_file is not None:
            raise ValueError("PMWave uses either seed or phase_file")
        if self.phase_generator not in ("numpy", "matlab"):
            raise ValueError("PMWave.phase_generator must be 'numpy' or 'matlab'")
        if self.phase_generator == "matlab" and self.phase_file is None and self.seed is None:
            raise ValueError("MATLAB phase generation needs a substream seed")
        if (not isinstance(self.discretization, str)
                or self.discretization.lower().replace("_", "")
                not in ("equalenergy", "traditional")):
            raise ValueError("PMWave.discretization must be 'equal_energy' or 'traditional'")
        if (self.directions is None) != (self.spreading is None):
            raise ValueError("PMWave directions and spreading must be supplied together")
        wave = {"type": "pm", "height": self.height,
                "period": self.period,
                "directions": (list(self.directions) if self.directions is not None
                               else [self.direction]),
                "spreading": (list(self.spreading) if self.spreading is not None
                              else [1.0]),
                "frequency_count": (self.frequency_count if self.frequency_count is not None
                                    else (1000 if self.discretization.lower() == "traditional"
                                          else 500))}
        if self.discretization.lower().replace("_", "") == "traditional":
            wave["discretization"] = "traditional"
        if self.frequency_range is not None:
            wave["frequency_range"] = list(self.frequency_range)
        if self.water_depth is not None:
            wave["water_depth"] = self.water_depth
        if self.current is not None:
            if not isinstance(self.current, Current):
                raise TypeError("PMWave.current must be Current")
            wave["current"] = self.current.as_case()
        if self.seed is not None:
            wave["seed"] = self.seed
        if self.phase_generator != "numpy":
            wave["phase_generator"] = self.phase_generator
        if self.phase_file is not None:
            wave["phase_file"] = str(self.phase_file)
        return wave


@dataclass(frozen=True)
class JONSWAPWave(PMWave):
    """JONSWAP sea using WEC-Sim's default or an explicit peak factor."""

    gamma: float | None = None

    def as_case(self) -> dict:
        wave = super().as_case()
        wave["type"] = "jonswap"
        if self.gamma is not None:
            wave["gamma"] = self.gamma
        return wave


@dataclass(frozen=True)
class MultiPMWave:
    """Sum independently phased PM seas in a fixed-hinge WEC.

    Each sea keeps its own height, period, heading, and phase seed or file.
    The published OSWEC multiple-spectrum case uses spline interpolation.
    """

    seas: tuple[PMWave, ...]
    excitation_interpolation: str = "linear"

    def as_case(self) -> dict:
        if len(self.seas) < 2:
            raise ValueError("MultiPMWave needs at least two PM seas")
        if self.excitation_interpolation not in ("linear", "spline"):
            raise ValueError("MultiPMWave excitation interpolation must be linear or spline")
        settings = []
        for sea in self.seas:
            if type(sea) is not PMWave:
                raise TypeError("MultiPMWave seas must be PMWave instances")
            setting = sea.as_case()
            if set(setting) - {"type", "height", "period", "directions",
                                "spreading", "seed", "phase_file",
                                "phase_generator", "frequency_count"}:
                raise ValueError("MultiPMWave supports equal-energy PM seas without current or depth overrides")
            setting.pop("type")
            settings.append(setting)
        return {"type": "pm_multi", "seas": settings,
                "excitation_interpolation": self.excitation_interpolation}


@dataclass(frozen=True)
class ImportedSpectrumWave:
    """WEC-Sim three-column MAT spectrum with its saved phase realization.

    ``file`` resolves from ``WEC.run(base_dir=...)``. The MAT file contains
    frequency in Hz, density in m²/Hz, and phase in radians.
    """

    file: str | Path

    def as_case(self) -> dict:
        return {"type": "spectrumImport", "file": str(self.file)}


@dataclass(frozen=True)
class FullDirectionalSpectrumWave:
    """Frequency-resolved directional MAT sea for a supported fixed hinge.

    Heading-bin width is included in excitation by default. The pinned
    MATLAB Full_Directional_Waves force block omitted it; select
    ``force_quadrature="matlab_omitted"`` only for that source comparison.
    Files resolve from ``WEC.run(base_dir=...)``.
    """

    file: str | Path
    seed: int | None = None
    phase_file: str | Path | None = None
    phase_generator: str = "numpy"
    excitation_interpolation: str = "linear"
    force_quadrature: str = "integrated"

    def as_case(self) -> dict:
        if not str(self.file):
            raise ValueError("full-directional wave needs a MAT file")
        if self.seed is not None and self.phase_file is not None:
            raise ValueError("full-directional wave uses either seed or phase_file")
        if self.phase_file is not None and not str(self.phase_file):
            raise ValueError("phase_file must be a nonempty path")
        if self.phase_generator not in ("numpy", "matlab"):
            raise ValueError("phase_generator must be 'numpy' or 'matlab'")
        if (self.phase_generator == "matlab" and self.phase_file is None
                and self.seed is None):
            raise ValueError("MATLAB phase generation needs a substream seed")
        if self.excitation_interpolation not in ("linear", "spline_frequency"):
            raise ValueError("excitation_interpolation must be linear or spline_frequency")
        if self.force_quadrature not in ("integrated", "matlab_omitted"):
            raise ValueError("force_quadrature must be integrated or matlab_omitted")
        wave = {"type": "spectrumImportFullDir", "file": str(self.file),
                "excitation_interpolation": self.excitation_interpolation,
                "force_quadrature": self.force_quadrature}
        if self.seed is not None:
            wave["seed"] = self.seed
        if self.phase_file is not None:
            wave["phase_file"] = str(self.phase_file)
        if self.phase_generator != "numpy":
            wave["phase_generator"] = self.phase_generator
        return wave


@dataclass(frozen=True)
class ImportedElevationWave:
    """Sampled time/elevation MAT record in WEC-Sim's ``elevationImport`` mode.

    ``file`` resolves from ``WEC.run(base_dir=...)``. The named MAT variable
    contains increasing time and elevation columns in seconds and metres.
    """

    file: str | Path
    variable: str = "etaData"
    reapply_force_ramp: bool = False

    def as_case(self) -> dict:
        return {"type": "elevationImport", "file": str(self.file),
                "variable": self.variable,
                "reapply_force_ramp": self.reapply_force_ramp}


@dataclass(frozen=True)
class NoWave:
    def as_case(self) -> dict:
        return {"type": "none"}

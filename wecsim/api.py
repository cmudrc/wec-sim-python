"""Python interface for configuring supported WEC devices.

The objects in this module assemble the same validated case used by the JSON
runner. Users can define bodies, named motions, attachment points, PTOs, and
waves directly in Python and receive NumPy arrays without writing JSON or CSV.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from numbers import Real
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np

from .caseDynamics import CaseResponse, run_case
from .controls import DeclutchingControl, LatchingControl
from .directLinearGenerator import DirectLinearGenerator
from .floatingOwc import FloatingOwcChamber, FloatingOwcTurbine
from .hardStops import LinearHardStops
from .morison import MorisonElement, finite_depth_wavenumber
from .moorDyn import MoorDyn
from .orifice import OrificePTO


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


@dataclass(frozen=True)
class MotionHistory:
    position: np.ndarray
    velocity: np.ndarray


@dataclass(frozen=True)
class FlexibleModeHistory:
    position: np.ndarray
    velocity: np.ndarray
    acceleration: np.ndarray


@dataclass(frozen=True)
class DirectDriveHistory:
    shaft_velocity: np.ndarray
    shaft_torque: np.ndarray
    inertia_torque: np.ndarray
    friction_torque: np.ndarray
    generator_torque: np.ndarray
    current: np.ndarray
    voltage: np.ndarray
    resistance_loss: np.ndarray
    electrical_power: np.ndarray
    mechanical_power: np.ndarray


@dataclass(frozen=True)
class LinearGeneratorHistory:
    flux_d: np.ndarray
    flux_q: np.ndarray
    electrical_angle: np.ndarray
    friction_force: np.ndarray
    electrical_power: np.ndarray
    phase_current: np.ndarray
    phase_voltage: np.ndarray


@dataclass(frozen=True)
class PTOHistory:
    """PTO stroke in metres, or angle in radians for a rotational PTO."""

    stroke: np.ndarray
    velocity: np.ndarray
    force: np.ndarray
    absorbed_power: np.ndarray
    direct_drive: DirectDriveHistory | None = None
    linear_generator: LinearGeneratorHistory | None = None


@dataclass(frozen=True)
class WECResult:
    time: np.ndarray
    bodies: dict[str, MotionHistory]
    coordinates: dict[str, MotionHistory]
    ptos: dict[str, PTOHistory]
    wave_elevation: np.ndarray | None
    case: dict
    raw: CaseResponse
    flexible_modes: dict[str, FlexibleModeHistory] = field(default_factory=dict)
    body_forces: dict[str, np.ndarray] = field(default_factory=dict)


class WEC:
    """Build supported WEC motion/PTO layouts or fixed Morison bodies.

    Mapped-coordinate rotations and PTO stroke are linearized about the
    reference pose. The fixed-hinge and floating-joint layouts use their
    dedicated kinematics. Hydrodynamic body order must match the body order
    in the supplied HDF5 file.
    """

    def __init__(self, name: str = "WEC", *, body_to_body: bool = False):
        if not isinstance(name, str) or not name:
            raise ValueError("WEC name must be a nonempty string")
        if not isinstance(body_to_body, bool):
            raise ValueError("body_to_body must be a boolean")
        self.name = name
        self.body_to_body = body_to_body
        self.bodies: list[Body] = []
        self.coordinates: list[Coordinate] = []
        self.ptos: list[LinearPTO] = []
        self.rotational_ptos: list[RotationalPTO] = []
        self._floating_gbm_body: Body | None = None
        self._floating_gbm_orifice: dict | None = None
        self._floating_joint: _FloatingJoint | None = None
        self._fixed_hinge: _FixedHinge | None = None
        self._floating_owc: _FloatingOwc | None = None
        self._morison_elements: list[tuple[Body, MorisonElement]] = []

    def body(self, name: str, hydro_file: str | Path, *,
             mass: str | float = "equilibrium",
             inertia: Sequence[float] = (0, 0, 0),
             hydro_body: int | None = None,
             mean_drift: str = "none",
             passive_yaw: bool = False,
             passive_yaw_threshold: float = 0.0,
             yaw_heading_bank: Sequence[float] | None = None,
             geometry_file: str | Path | None = None,
             nonlinear_hydro: str | None = None,
             drag_coefficient: float = 0.0,
             drag_area: float = 0.0) -> Body:
        if any(existing.name == name for existing in self.bodies):
            raise ValueError(f"body name already exists: {name}")
        if (isinstance(passive_yaw_threshold, bool)
                or not isinstance(passive_yaw_threshold, Real)
                or not np.isfinite(passive_yaw_threshold)
                or passive_yaw_threshold < 0
                or (passive_yaw_threshold > 0 and passive_yaw is not True)):
            raise ValueError("passive_yaw_threshold needs a nonnegative degree value and passive_yaw=True")
        if yaw_heading_bank is not None:
            headings = np.asarray(yaw_heading_bank, dtype=float)
            if (not passive_yaw or passive_yaw_threshold
                    or headings.ndim != 1 or len(headings) < 2
                    or not np.isfinite(headings).all()
                    or not np.all(np.diff(headings) > 0)
                    or headings[0] < -180 or headings[-1] > 180):
                raise ValueError("yaw_heading_bank needs passive_yaw, no threshold, and ordered directions in [-180, 180]")
            yaw_heading_bank = tuple(float(value) for value in headings)
        body = Body(name, hydro_file, mass, tuple(inertia), hydro_body,
                    mean_drift, passive_yaw, geometry_file, nonlinear_hydro,
                    drag_coefficient, drag_area,
                    passive_yaw_threshold=passive_yaw_threshold,
                    yaw_heading_bank=yaw_heading_bank)
        self.bodies.append(body)
        return body

    def variable_body(self, name: str, states: Sequence[HydroState], *,
                      switch_times: Sequence[float]) -> Body:
        """Add a body whose draft, mass, and BEM data switch together.

        Current dynamics support one heave-only body in regular waves.
        ``switch_times`` must contain one grid-aligned time per transition.
        """
        ordered = tuple(states)
        if (len(ordered) < 2 or not all(isinstance(s, HydroState) for s in ordered)
                or len(switch_times) != len(ordered) - 1):
            raise ValueError("variable body needs states and one fewer switch times")
        if any(existing.name == name for existing in self.bodies):
            raise ValueError(f"body name already exists: {name}")
        first = ordered[0]
        body = Body(name, first.hydro_file, first.mass, first.inertia,
                    variable_hydro=VariableHydro(ordered, tuple(switch_times)))
        self.bodies.append(body)
        return body

    def fixed_body(self, name: str, *, center_gravity: Sequence[float],
                   mass: str | float = "equilibrium",
                   inertia: Sequence[float] = (0, 0, 0),
                   volume: float = 0.0) -> Body:
        """Add a stationary body without HDF5 hydrodynamics."""
        if any(existing.name == name for existing in self.bodies):
            raise ValueError(f"body name already exists: {name}")
        center = tuple(center_gravity)
        if len(center) != 3 or not np.isfinite(center).all():
            raise ValueError("center_gravity needs three finite coordinates")
        body = Body(name, None, mass, tuple(inertia), fixed=True,
                    center_gravity=center, volume=volume)
        self.bodies.append(body)
        return body

    def morison_element(self, body: Body, *, point: BodyPoint,
                        drag_coefficient: Sequence[float],
                        added_mass_coefficient: Sequence[float],
                        area: Sequence[float], volume: float,
                        phase_mode: str = "directional") -> MorisonElement:
        """Attach a Cartesian Morison element at a body-local point.

        Moving hydrodynamic bodies support axial elements in pure heave or
        regular-wave surge/heave/pitch motion; the case runner validates
        those layouts.
        """
        if not any(body is item for item in self.bodies):
            raise ValueError("Morison body must belong to this WEC")
        if not isinstance(point, BodyPoint) or point.body is not body:
            raise ValueError("Morison point must belong to its body")
        element = MorisonElement(
            tuple(point.coordinates()), tuple(drag_coefficient),
            tuple(added_mass_coefficient), tuple(area), volume, phase_mode,
        )
        self._morison_elements.append((body, element))
        return element

    def floating_gbm(
        self, body: Body, *, orifice: OrificePTO | None = None,
        pto_name: str = "orifice",
        orifice_force_path: str = "coupled",
        heave_linear_damping: float = 0,
        mode_linear_damping: float = 0,
        heave_drag_cd: float = 0, heave_drag_area: float = 0,
        pitch_drag_cd: float = 0, pitch_drag_area: float = 0,
    ) -> None:
        """Select a floating surge/heave/pitch joint with HDF5 flexible modes.

        The joint and body center of gravity must coincide at the origin.
        Without an orifice, support is a regular-wave body without a PTO.
        With an orifice, the one-mode OWC layout supports a PM sea. The default
        couples the piston to rigid heave and the flexible state. The optional
        ``published_owc`` path reproduces the source's omitted flexible-state
        piston force for paired comparisons.
        """
        if not any(body is item for item in self.bodies):
            raise ValueError("floating GBM body must belong to this WEC")
        if self._floating_gbm_body is not None or self._floating_owc is not None:
            raise ValueError("a floating GBM body is already selected")
        if orifice is not None and not isinstance(orifice, OrificePTO):
            raise TypeError("orifice must be an OrificePTO")
        if not isinstance(orifice_force_path, str) or orifice_force_path not in (
                "coupled", "published_owc"):
            raise ValueError("orifice_force_path must be coupled or published_owc")
        if orifice is None and orifice_force_path != "coupled":
            raise ValueError("published_owc force path needs an orifice")
        settings = dict(
            heave_linear_damping=heave_linear_damping,
            mode_linear_damping=mode_linear_damping,
            heave_drag_cd=heave_drag_cd,
            heave_drag_area=heave_drag_area,
            pitch_drag_cd=pitch_drag_cd,
            pitch_drag_area=pitch_drag_area,
        )
        if (any(isinstance(value, bool) or not isinstance(value, Real)
                for value in settings.values())
                or not np.isfinite(list(settings.values())).all()
                or any(value < 0 for value in settings.values())):
            raise ValueError("floating GBM damping and drag settings must be nonnegative and finite")
        if orifice is None and any(settings.values()):
            raise ValueError("floating GBM damping and drag currently need an orifice")
        if not isinstance(pto_name, str) or not pto_name:
            raise ValueError("pto_name must be a nonempty string")
        self._floating_gbm_body = body
        if orifice is not None:
            self._floating_gbm_orifice = {"orifice": orifice, "name": pto_name,
                                          "orifice_force_path": orifice_force_path,
                                          **settings}

    def floating_joint(self, float_body: Body, spar_body: Body, *,
                       location: WorldPoint = WorldPoint(0, 0, 0),
                       pto_name: str = "relative_heave",
                       damping: float = 0.0, stiffness: float = 0.0,
                       equilibrium_position: float = 0.0,
                       mooring_surge_stiffness: float = 0.0,
                       moordyn: MoorDyn | None = None,
                       moordyn_point: BodyPoint | None = None,
                       hard_stops: LinearHardStops | None = None,
                       radiation_method: str | None = None,
                       added_mass_scheme: str = "implicit") -> None:
        """Select the paired two-body surge/heave/pitch slider joint.

        The PTO acts on float heave minus spar heave. The reduced joint does
        not model off-axis PTO endpoints or arbitrary Simscape constraints.
        A native MoorDyn model can attach to one spar-local ``BodyPoint``;
        it is advanced once per time step with convolution radiation.
        """
        if (len(self.bodies) != 2 or self.bodies[0] is not float_body
                or self.bodies[1] is not spar_body or self._floating_joint is not None
                or self._floating_gbm_body is not None or self._fixed_hinge is not None
                or self._floating_owc is not None
                or self.coordinates
                or self.ptos or self.rotational_ptos or self._morison_elements):
            raise ValueError("floating_joint needs exactly two ordered bodies and no other layout or PTO")
        if not isinstance(location, WorldPoint) or location.x != 0 or location.y != 0:
            raise ValueError("floating_joint location needs a world point on the z axis")
        if not isinstance(pto_name, str) or not pto_name:
            raise ValueError("floating_joint PTO name must be nonempty")
        if hard_stops is not None and not isinstance(hard_stops, LinearHardStops):
            raise TypeError("hard_stops must be LinearHardStops")
        if radiation_method not in (None, "constant", "convolution", "fir"):
            raise ValueError("unsupported floating_joint radiation method")
        if added_mass_scheme not in ("implicit", "simulink_delay"):
            raise ValueError("unsupported floating_joint added-mass scheme")
        if moordyn is not None:
            if (not isinstance(moordyn, MoorDyn)
                    or not isinstance(moordyn_point, BodyPoint)
                    or moordyn_point.body is not spar_body
                    or not np.isfinite(moordyn_point.coordinates()).all()):
                raise ValueError("MoorDyn needs a session and a finite spar-local point")
            if mooring_surge_stiffness:
                raise ValueError("MoorDyn and the joint surge spring cannot be combined")
        elif moordyn_point is not None:
            raise ValueError("a MoorDyn attachment point needs a MoorDyn session")
        self._floating_joint = _FloatingJoint(
            float_body, spar_body, location, pto_name, damping, stiffness,
            equilibrium_position, mooring_surge_stiffness, hard_stops,
            radiation_method, added_mass_scheme, moordyn, moordyn_point,
        )

    def floating_owc(
        self, floater: Body, column: Body, *, moordyn: MoorDyn,
        column_height: float, column_diameter: float,
        moordyn_point: BodyPoint | None = None,
        pto_name: str = "column_slider", pto_stiffness: float = 0.0,
        pto_damping: float = 0.0,
        chamber: FloatingOwcChamber | None = None,
        turbine: FloatingOwcTurbine | None = None,
        initial_rotor_speed: float = 150.0,
    ) -> None:
        """Select the coupled two-body OWC, axial PTO, air train, and mooring.

        The floater has six motions and the water column slides on its local
        z axis. Both bodies must use the same two-body HDF5 file in order.
        """
        if (len(self.bodies) != 2 or self.bodies[0] is not floater
                or self.bodies[1] is not column or self.body_to_body
                or self._floating_owc is not None or self._floating_joint is not None
                or self._fixed_hinge is not None or self._floating_gbm_body is not None
                or self.coordinates or self.ptos or self.rotational_ptos
                or self._morison_elements):
            raise ValueError("floating_owc needs exactly two ordered bodies and no other layout")
        if not isinstance(moordyn, MoorDyn):
            raise TypeError("floating_owc needs a MoorDyn session")
        if moordyn_point is not None and (
                not isinstance(moordyn_point, BodyPoint)
                or moordyn_point.body is not floater
                or not np.isfinite(moordyn_point.coordinates()).all()):
            raise ValueError("MoorDyn point must be finite and floater-local")
        if chamber is not None and not isinstance(chamber, FloatingOwcChamber):
            raise TypeError("chamber must be FloatingOwcChamber")
        if turbine is not None and not isinstance(turbine, FloatingOwcTurbine):
            raise TypeError("turbine must be FloatingOwcTurbine")
        if not isinstance(pto_name, str) or not pto_name:
            raise ValueError("floating_owc PTO name must be nonempty")
        settings = (column_height, column_diameter, pto_stiffness,
                    pto_damping, initial_rotor_speed)
        if (any(isinstance(value, bool) or not isinstance(value, Real)
                for value in settings) or not np.isfinite(settings).all()
                or column_height <= 0 or column_diameter <= 0
                or initial_rotor_speed <= 0):
            raise ValueError("floating_owc settings must be finite and physical")
        self._floating_owc = _FloatingOwc(
            floater, column, moordyn, moordyn_point,
            float(column_height), float(column_diameter), pto_name,
            float(pto_stiffness), float(pto_damping), chamber, turbine,
            float(initial_rotor_speed),
        )

    def fixed_hinge(self, flap: Body, base: Body | None = None, *,
                    location: WorldPoint = WorldPoint(0, 0, -10),
                    pto_location: WorldPoint = WorldPoint(0, 0, -8.9),
                    pto_name: str = "hinge", damping: float = 0.0,
                    stiffness: float = 0.0, equilibrium_angle: float = 0.0,
                    added_mass_scheme: str = "implicit") -> None:
        """Select the published OSWEC pitch hinge and torsional PTO.

        The flap moves in surge, heave, and pitch about a fixed base. Both
        joint and PTO points currently lie on the world z axis.
        """
        if (len(self.bodies) != (1 if base is None else 2)
                or self.bodies[0] is not flap
                or base is not None and self.bodies[1] is not base
                or flap.fixed or self.body_to_body
                or self._fixed_hinge is not None or self._floating_joint is not None
                or self._floating_owc is not None
                or self._floating_gbm_body is not None or self.coordinates
                or self.ptos or self.rotational_ptos or self._morison_elements):
            raise ValueError("fixed_hinge needs one flap, optional ordered base, and no other layout or PTO")
        if (not isinstance(location, WorldPoint)
                or not isinstance(pto_location, WorldPoint)
                or not np.isfinite(location.coordinates()).all()
                or not np.isfinite(pto_location.coordinates()).all()
                or location.x != 0 or location.y != 0
                or pto_location.x != 0 or pto_location.y != 0):
            raise ValueError("fixed_hinge points need finite world locations on the z axis")
        if not isinstance(pto_name, str) or not pto_name:
            raise ValueError("fixed_hinge PTO name must be nonempty")
        if (not np.isfinite([damping, stiffness, equilibrium_angle]).all()
                or damping < 0 or stiffness < 0
                or equilibrium_angle and not stiffness):
            raise ValueError("fixed_hinge PTO needs finite nonnegative damping and stiffness")
        if added_mass_scheme not in ("implicit", "simulink_delay"):
            raise ValueError("unsupported fixed_hinge added-mass scheme")
        self._fixed_hinge = _FixedHinge(
            flap, base, location, pto_location, pto_name,
            damping, stiffness, equilibrium_angle, added_mass_scheme,
        )

    def coordinate(self, name: str, *motions: Motion) -> Coordinate:
        if any(existing.name == name for existing in self.coordinates):
            raise ValueError(f"coordinate name already exists: {name}")
        if not motions:
            raise ValueError("a coordinate needs body motions")
        if any(not isinstance(motion, Motion)
               or not any(motion.body is body for body in self.bodies)
               for motion in motions):
            raise ValueError("coordinate motions must use bodies in this WEC")
        coordinate = Coordinate(name, tuple(motions))
        self.coordinates.append(coordinate)
        return coordinate

    def pto(self, name: str, from_point: BodyPoint | WorldPoint,
            to_point: BodyPoint | WorldPoint, *,
            axis: Sequence[float] | None = None,
            damping: float = 0.0, stiffness: float = 0.0,
            equilibrium_position: float | None = None,
            pretension: float | None = None,
            control: DeclutchingControl | LatchingControl | None = None,
            direct_drive: SimpleDirectDrive | None = None,
            linear_generator: DirectLinearGenerator | None = None) -> LinearPTO:
        if any(existing.name == name for existing in
               (*self.ptos, *self.rotational_ptos)):
            raise ValueError(f"PTO name already exists: {name}")
        for point in (from_point, to_point):
            if not isinstance(point, (BodyPoint, WorldPoint)):
                raise TypeError("PTO endpoints must be body points or world points")
            if (isinstance(point, BodyPoint)
                    and not any(point.body is body for body in self.bodies)):
                raise ValueError("PTO body point must belong to this WEC")
        if control is not None:
            if not isinstance(control, (DeclutchingControl, LatchingControl)):
                raise TypeError("PTO control must be a supported control law")
            if damping or stiffness or equilibrium_position is not None or pretension is not None:
                raise ValueError("sampled PTO control sets its own gain and has no spring")
        if direct_drive is not None:
            if not isinstance(direct_drive, SimpleDirectDrive):
                raise TypeError("direct_drive must be SimpleDirectDrive")
            if (control is not None or damping or stiffness
                    or equilibrium_position is not None or pretension is not None
                    or linear_generator is not None):
                raise ValueError("direct drive supplies its own PTO force")
        if linear_generator is not None:
            if not isinstance(linear_generator, DirectLinearGenerator):
                raise TypeError("linear_generator must be DirectLinearGenerator")
            if (control is not None or damping or stiffness
                    or equilibrium_position is not None or pretension is not None):
                raise ValueError("linear generator supplies its own PTO force")
        pto = LinearPTO(
            name, from_point, to_point,
            tuple(axis) if axis is not None else None,
            damping, stiffness, equilibrium_position, pretension, control,
            direct_drive, linear_generator,
        )
        self.ptos.append(pto)
        return pto

    def rotational_pto(self, name: str, coordinate: Coordinate, *,
                       damping: float = 0.0, stiffness: float = 0.0,
                       equilibrium_angle: float = 0.0) -> RotationalPTO:
        """Attach a torsional PTO to a named rotation coordinate.

        Damping is in N m s/rad, stiffness in N m/rad, and equilibrium angle
        in radians. The returned PTO history reports angle as ``stroke``.
        """
        if not any(coordinate is item for item in self.coordinates):
            raise ValueError("rotational PTO coordinate must belong to this WEC")
        if any(existing.name == name for existing in
               (*self.ptos, *self.rotational_ptos)):
            raise ValueError(f"PTO name already exists: {name}")
        pto = RotationalPTO(name, coordinate, damping, stiffness,
                            equilibrium_angle)
        self.rotational_ptos.append(pto)
        return pto

    def to_case(
        self, wave: RegularWave | RegularCICWave | PMWave | JONSWAPWave | ImportedSpectrumWave | FullDirectionalSpectrumWave | ImportedElevationWave | NoWave, *,
        dt: float, end_time: float,
        ramp_time: float | None = None,
        radiation_memory: float | None = None,
        rho: float | None = None, g: float | None = None,
        initial_coordinate: Mapping[str, float] | Sequence[float] | None = None,
        initial_speed: Mapping[str, float] | Sequence[float] | None = None,
    ) -> dict:
        """Return the case mapping used by the validated dynamics runner."""
        if not isinstance(wave, (RegularWave, RegularCICWave, PMWave,
                                 JONSWAPWave, ImportedSpectrumWave,
                                 FullDirectionalSpectrumWave,
                                 ImportedElevationWave, NoWave)):
            raise TypeError("wave must be a supported regular, irregular, imported, or no-wave configuration")
        simulation = {"dt": dt, "end_time": end_time}
        for key, value in (
            ("ramp_time", ramp_time), ("radiation_memory", radiation_memory),
            ("rho", rho), ("g", g),
        ):
            if value is not None:
                simulation[key] = value
        if self._floating_joint is not None:
            return self._floating_joint_case(
                wave, simulation, initial_coordinate, initial_speed,
            )
        if self._floating_owc is not None:
            return self._floating_owc_case(
                wave, simulation, radiation_memory,
                initial_coordinate, initial_speed,
            )
        if self._fixed_hinge is not None:
            return self._fixed_hinge_case(
                wave, simulation, initial_coordinate, initial_speed,
            )
        bodies = []
        for index, body in enumerate(self.bodies, start=1):
            if body.fixed:
                body_case = {
                    "name": body.name, "nonhydro": True, "fixed": True,
                    "center_gravity": list(body.center_gravity),
                    "mass": body.mass, "inertia": list(body.inertia),
                    "volume": body.volume,
                    "morison_elements": [
                        {"point": list(element.point),
                         "drag_coefficient": list(element.drag_coefficient),
                         "added_mass_coefficient": list(element.added_mass_coefficient),
                         "area": list(element.area), "volume": element.volume,
                         "phase_mode": element.phase_mode}
                        for attached_body, element in self._morison_elements
                        if attached_body is body
                    ],
                }
                bodies.append(body_case)
                continue
            body_case = {
                "name": body.name,
                "hydro_file": str(body.hydro_file),
                "hydro_body": body.hydro_body if body.hydro_body is not None else index,
                "mass": body.mass,
                "inertia": list(body.inertia),
            }
            if body.mean_drift != "none":
                body_case["mean_drift"] = body.mean_drift
            if body.passive_yaw:
                body_case["passive_yaw"] = True
            if body.passive_yaw_threshold:
                body_case["passive_yaw_threshold"] = body.passive_yaw_threshold
            if body.yaw_heading_bank is not None:
                body_case["yaw_heading_bank"] = list(body.yaw_heading_bank)
            if body.geometry_file is not None:
                body_case["geometry_file"] = str(body.geometry_file)
            if body.nonlinear_hydro is not None:
                body_case["nonlinear_hydro"] = body.nonlinear_hydro
            if body.drag_coefficient or body.drag_area:
                body_case["drag_coefficient"] = body.drag_coefficient
                body_case["drag_area"] = body.drag_area
            if body.variable_hydro is not None:
                body_case["variable_hydro"] = {
                    "states": [
                        {"hydro_file": str(state.hydro_file),
                         "mass": state.mass,
                         "inertia": list(state.inertia)}
                        for state in body.variable_hydro.states
                    ],
                    "switch_times": list(body.variable_hydro.switch_times),
                }
            elements = [element for attached, element in self._morison_elements
                        if attached is body]
            if elements:
                body_case["morison_elements"] = [
                    {"point": list(element.point),
                     "drag_coefficient": list(element.drag_coefficient),
                     "added_mass_coefficient": list(element.added_mass_coefficient),
                     "area": list(element.area), "volume": element.volume,
                     "phase_mode": element.phase_mode}
                    for element in elements
                ]
            bodies.append(body_case)
        constraint = {"kind": "linear_subspace", "coordinates": []}
        if any(body.fixed for body in self.bodies):
            if (not all(body.fixed for body in self.bodies)
                    or self.coordinates or self.ptos or self.rotational_ptos
                    or self.body_to_body or radiation_memory is not None
                    or initial_coordinate is not None or initial_speed is not None):
                raise ValueError("fixed Morison bodies need no moving coordinates, PTOs, or radiation")
            constraint = {"kind": "fixed_morison"}
        if self._floating_gbm_body is not None:
            if (len(self.bodies) != 1 or self.bodies[0] is not self._floating_gbm_body
                    or self.coordinates or self.ptos or self.rotational_ptos
                    or self.body_to_body
                    or initial_coordinate is not None or initial_speed is not None
                    or (self._floating_gbm_orifice is None
                        and (not isinstance(wave, RegularWave)
                             or radiation_memory is not None))
                    or (self._floating_gbm_orifice is not None
                        and type(wave) is not PMWave)):
                raise ValueError("floating GBM needs one body and a supported regular or PM/orifice layout")
            constraint = {"kind": "floating_gbm", "location": [0, 0, 0]}
            if self._floating_gbm_orifice is not None:
                config = self._floating_gbm_orifice
                constraint["orifice"] = {"name": config["name"],
                                          **asdict(config["orifice"])}
                constraint.update({key: value for key, value in config.items()
                                   if key not in ("name", "orifice")})
        for coordinate in self.coordinates:
            motions = []
            for motion in coordinate.motions:
                entry = {"body": motion.body.name, "dof": motion.dof,
                         "scale": motion.scale}
                if isinstance(motion.pivot, WorldPoint):
                    entry["pivot"] = {"world": motion.pivot.coordinates()}
                elif isinstance(motion.pivot, BodyPoint):
                    if motion.pivot.body is not motion.body:
                        raise ValueError("body-local pivot must belong to its motion body")
                    entry["pivot"] = {"point": motion.pivot.coordinates()}
                elif motion.pivot is not None:
                    raise TypeError("pivot must be a world or body point")
                motions.append(entry)
            constraint["coordinates"].append({"name": coordinate.name,
                                               "motions": motions})
        if initial_coordinate is not None:
            constraint["initial_coordinate"] = _state(initial_coordinate)
        if initial_speed is not None:
            constraint["initial_speed"] = _state(initial_speed)
        case = {
            "name": self.name,
            "simulation": simulation,
            "wave": wave.as_case(),
            "bodies": bodies,
            "constraint": constraint,
        }
        if self.body_to_body:
            case["body_to_body"] = True
        if self.ptos or self.rotational_ptos:
            case["ptos"] = ([self._pto_case(pto) for pto in self.ptos]
                            + [self._rotational_pto_case(pto)
                               for pto in self.rotational_ptos])
        return case

    def _floating_joint_case(self, wave, simulation,
                             initial_coordinate, initial_speed) -> dict:
        joint = self._floating_joint
        if (len(self.bodies) != 2 or self.bodies[0] is not joint.float_body
                or self.bodies[1] is not joint.spar_body or self.coordinates
                or self.ptos or self.rotational_ptos or self._morison_elements
                or self._floating_gbm_body is not None):
            raise ValueError("floating_joint cannot combine with other bodies, coordinates, or PTOs")
        if not isinstance(wave, (RegularWave, RegularCICWave, PMWave,
                                 ImportedElevationWave, NoWave)):
            raise ValueError("floating_joint supports regular, irregular, imported elevation, or no waves")
        if (not np.isfinite(joint.location.coordinates()).all()
                or not np.isfinite([joint.damping, joint.stiffness,
                                    joint.equilibrium_position,
                                    joint.mooring_surge_stiffness]).all()
                or joint.mooring_surge_stiffness < 0):
            raise ValueError("floating_joint location, PTO, and mooring settings must be finite")
        bodies = []
        for index, body in enumerate((joint.float_body, joint.spar_body), start=1):
            if (body.fixed or body.hydro_file is None or body.mass != "equilibrium"
                    or body.mean_drift != "none" or body.passive_yaw
                    or body.geometry_file is not None or body.nonlinear_hydro is not None
                    or body.variable_hydro is not None or body.drag_coefficient
                    or body.drag_area or body.passive_yaw_threshold
                    or body.yaw_heading_bank is not None):
                raise ValueError("floating_joint needs equilibrium-mass hydrodynamic bodies without extra force models")
            inertia = np.asarray(body.inertia, dtype=float)
            if inertia.shape != (3,) or not np.isfinite(inertia).all() or inertia[1] <= 0:
                raise ValueError("floating_joint needs a positive pitch inertia on each body")
            bodies.append({
                "hydro_file": str(body.hydro_file),
                "hydro_body": body.hydro_body if body.hydro_body is not None else index,
                "mass": "equilibrium", "pitch_inertia": float(inertia[1]),
            })
        simulation["added_mass_scheme"] = joint.added_mass_scheme
        if joint.radiation_method is not None:
            simulation["radiation_method"] = joint.radiation_method
        constraint = {"kind": "floating_joint",
                      "location": joint.location.coordinates()}
        if initial_coordinate is not None:
            constraint["initial_coordinate"] = _state(initial_coordinate)
        if initial_speed is not None:
            constraint["initial_speed"] = _state(initial_speed)
        pto = {"kind": "relative_heave", "damping": joint.damping,
               "stiffness": joint.stiffness,
               "equilibrium_position": joint.equilibrium_position}
        if joint.hard_stops is not None:
            pto["hard_stops"] = asdict(joint.hard_stops)
        case = {"name": self.name, "simulation": simulation,
                "wave": wave.as_case(), "bodies": bodies,
                "constraint": constraint, "pto": pto}
        if self.body_to_body:
            case["body_to_body"] = True
        if joint.mooring_surge_stiffness:
            case["mooring"] = {"kind": "joint_surge_spring",
                               "stiffness": joint.mooring_surge_stiffness}
        if joint.moordyn is not None:
            case["mooring"] = {"kind": "moor_dyn", "session": joint.moordyn,
                               "point": joint.moordyn_point.coordinates()}
        return case

    def _floating_owc_case(self, wave, simulation, radiation_memory,
                           initial_coordinate, initial_speed) -> dict:
        layout = self._floating_owc
        if (len(self.bodies) != 2 or self.bodies[0] is not layout.floater
                or self.bodies[1] is not layout.column
                or self._floating_joint is not None or self._fixed_hinge is not None
                or self._floating_gbm_body is not None or self.coordinates
                or self.ptos or self.rotational_ptos or self._morison_elements
                or type(wave) is not RegularWave or wave.direction != 0
                or radiation_memory is not None or initial_coordinate is not None
                or initial_speed is not None):
            raise ValueError("floating_owc currently needs zero-heading regular waves and no radiation memory or custom initial state")
        if (layout.floater.hydro_body not in (None, 1)
                or layout.column.hydro_body not in (None, 2)):
            raise ValueError("floating_owc bodies need ordered HDF5 body numbers")
        for body in (layout.floater, layout.column):
            if (body.fixed or body.hydro_file is None or body.mean_drift != "none"
                    or body.passive_yaw or body.geometry_file is not None
                    or body.nonlinear_hydro is not None
                    or body.variable_hydro is not None or body.drag_coefficient
                    or body.drag_area or body.passive_yaw_threshold
                    or body.yaw_heading_bank is not None):
                raise ValueError("floating_owc bodies cannot add other force models")
        if layout.floater.mass != "equilibrium":
            if (isinstance(layout.floater.mass, bool)
                    or not isinstance(layout.floater.mass, Real)
                    or not np.isfinite(layout.floater.mass)
                    or layout.floater.mass <= 0):
                raise ValueError("floater mass needs equilibrium or a positive value")
        if (isinstance(layout.column.mass, bool)
                or not isinstance(layout.column.mass, Real)
                or not np.isfinite(layout.column.mass)
                or layout.column.mass <= 0):
            raise ValueError("column mass must be positive")
        floater_inertia = np.asarray(layout.floater.inertia, dtype=float)
        column_inertia = np.asarray(layout.column.inertia, dtype=float)
        if (floater_inertia.shape != (3,) or column_inertia.shape != (3,)
                or not np.isfinite(floater_inertia).all()
                or not np.isfinite(column_inertia).all()
                or np.any(floater_inertia <= 0)
                or (np.any(column_inertia != 0) and np.any(column_inertia <= 0))):
            raise ValueError("floater inertia must be positive; column inertia is positive or automatic")
        simulation.setdefault("ramp_time", 50.0)
        return {
            "name": self.name, "simulation": simulation,
            "wave": wave.as_case(),
            "bodies": [
                {"name": layout.floater.name, "hydro_file": str(layout.floater.hydro_file),
                 "hydro_body": 1, "mass": layout.floater.mass,
                 "inertia": floater_inertia.tolist()},
                {"name": layout.column.name, "hydro_file": str(layout.column.hydro_file),
                 "hydro_body": 2, "mass": float(layout.column.mass),
                 "inertia": column_inertia.tolist()},
            ],
            "constraint": {
                "kind": "floating_owc", "column_height": layout.column_height,
                "column_diameter": layout.column_diameter,
                "initial_rotor_speed": layout.initial_rotor_speed,
                "chamber": layout.chamber, "turbine": layout.turbine,
            },
            "pto": {"kind": "column_slider", "name": layout.pto_name,
                    "stiffness": layout.pto_stiffness,
                    "damping": layout.pto_damping},
            "mooring": {"kind": "moor_dyn", "session": layout.moordyn,
                        "point": (layout.moordyn_point.coordinates()
                                  if layout.moordyn_point is not None else None)},
        }

    def _fixed_hinge_case(self, wave, simulation,
                          initial_coordinate, initial_speed) -> dict:
        hinge = self._fixed_hinge
        if (len(self.bodies) != (1 if hinge.base is None else 2)
                or self.bodies[0] is not hinge.flap
                or hinge.base is not None and self.bodies[1] is not hinge.base
                or self.coordinates or self.ptos or self.rotational_ptos
                or self._morison_elements or self._floating_gbm_body is not None
                or self._floating_joint is not None or self.body_to_body
                or initial_coordinate is not None or initial_speed is not None):
            raise ValueError("fixed_hinge needs its selected bodies and no custom coordinates or PTOs")
        if not isinstance(wave, (RegularWave, PMWave, FullDirectionalSpectrumWave)):
            raise ValueError("fixed_hinge supports regular, PM, or full-directional waves")
        if (isinstance(wave, RegularWave)
                and (wave.direction != 0
                     or hinge.added_mass_scheme != "implicit")):
            raise ValueError("regular fixed_hinge needs zero-heading waves and implicit added mass")
        if hinge.base is not None and hinge.base.fixed:
            if (hinge.base.center_gravity is None
                    or not np.isfinite(hinge.base.center_gravity).all()):
                raise ValueError("fixed nonhydrodynamic base needs a finite center of gravity")
            base_case = {
                "name": hinge.base.name, "nonhydro": True, "fixed": True,
                "center_gravity": list(hinge.base.center_gravity),
                "mass": hinge.base.mass, "inertia": list(hinge.base.inertia),
                "volume": hinge.base.volume,
            }
        elif hinge.base is not None:
            base_case = {"name": hinge.base.name,
                         "hydro_file": str(hinge.base.hydro_file),
                         "hydro_body": hinge.base.hydro_body or 2,
                         "mass": hinge.base.mass,
                         "pitch_inertia": float(hinge.base.inertia[1])}
        else:
            base_case = None
        flap = hinge.flap
        if (not isinstance(flap.mass, Real) or isinstance(flap.mass, bool)
                or not np.isfinite(flap.mass) or flap.mass <= 0
                or len(flap.inertia) != 3
                or not np.isfinite(flap.inertia).all()
                or flap.inertia[1] <= 0
                or flap.passive_yaw or flap.yaw_heading_bank is not None
                or flap.variable_hydro is not None
                or flap.nonlinear_hydro not in (None, "instantaneous")
                or (flap.geometry_file is None) != (flap.nonlinear_hydro is None)
                or flap.mean_drift != "none" or flap.drag_coefficient
                or flap.drag_area):
            raise ValueError("fixed_hinge needs positive flap mass and pitch inertia with supported mesh hydro settings")
        bodies = [{"name": flap.name, "hydro_file": str(flap.hydro_file),
                   "hydro_body": flap.hydro_body or 1,
                   "mass": float(flap.mass),
                   "pitch_inertia": float(flap.inertia[1])}]
        if flap.nonlinear_hydro is not None:
            bodies[0]["nonlinear_hydro"] = flap.nonlinear_hydro
            bodies[0]["geometry_file"] = str(flap.geometry_file)
        if base_case is not None:
            bodies.append(base_case)
        simulation["added_mass_scheme"] = hinge.added_mass_scheme
        return {
            "name": self.name, "simulation": simulation,
            "wave": wave.as_case(), "bodies": bodies,
            "constraint": {"kind": "fixed_hinge",
                           "location": hinge.location.coordinates()},
            "pto": {"kind": "pitch", "damping": hinge.damping,
                    "stiffness": hinge.stiffness,
                    "equilibrium_position": hinge.equilibrium_angle,
                    "location": hinge.pto_location.coordinates()},
        }

    def run(
        self, wave: RegularWave | RegularCICWave | PMWave | JONSWAPWave | ImportedSpectrumWave | FullDirectionalSpectrumWave | ImportedElevationWave | NoWave, *,
        dt: float, end_time: float,
        ramp_time: float | None = None,
        radiation_memory: float | None = None,
        rho: float | None = None, g: float | None = None,
        initial_coordinate: Mapping[str, float] | Sequence[float] | None = None,
        initial_speed: Mapping[str, float] | Sequence[float] | None = None,
        base_dir: str | Path = ".",
    ) -> WECResult:
        """Simulate the WEC and return named NumPy motion and PTO histories."""
        case = self.to_case(
            wave, dt=dt, end_time=end_time, ramp_time=ramp_time,
            radiation_memory=radiation_memory, rho=rho, g=g,
            initial_coordinate=initial_coordinate, initial_speed=initial_speed,
        )
        response = run_case(case, base_dir=base_dir)
        extras = dict(response.extra_outputs)
        bodies = {
            body.name: MotionHistory(
                response.body_position[:, index, :],
                response.body_velocity[:, index, :],
            )
            for index, body in enumerate(self.bodies)
        }
        if self._floating_joint is not None:
            coordinates = {
                name: MotionHistory(
                    response.coordinate_position[:, index],
                    response.coordinate_velocity[:, index],
                )
                for index, name in enumerate(
                    ("surge", "float_heave", "spar_heave", "pitch")
                )
            }
            pto = PTOHistory(
                response.pto_stroke,
                response.pto_velocity,
                response.pto_force,
                response.pto_absorbed_power,
            )
            return WECResult(
                response.time, bodies, coordinates,
                {self._floating_joint.pto_name: pto},
                response.wave_elevation, case, response,
            )
        if self._floating_owc is not None:
            pto = PTOHistory(
                response.pto_stroke, response.pto_velocity,
                response.pto_force, response.pto_absorbed_power,
            )
            return WECResult(
                response.time, bodies,
                {"column_stroke": MotionHistory(response.pto_stroke,
                                                response.pto_velocity)},
                {self._floating_owc.pto_name: pto},
                response.wave_elevation, case, response,
            )
        if self._fixed_hinge is not None:
            angle = response.body_position[:, 0, 4]
            speed = response.body_velocity[:, 0, 4]
            pto = PTOHistory(
                angle, speed, response.pto_force,
                self._fixed_hinge.damping * speed**2,
            )
            return WECResult(
                response.time, bodies,
                {"pitch": MotionHistory(angle, speed)},
                {self._fixed_hinge.pto_name: pto},
                response.wave_elevation, case, response,
            )
        coordinates = {
            coordinate.name: MotionHistory(
                extras[f"coordinate_{coordinate.name}_position"],
                extras[f"coordinate_{coordinate.name}_velocity"],
            )
            for coordinate in self.coordinates
        }
        ptos = {
            pto.name: PTOHistory(
                extras[f"pto_{pto.name}_stroke"],
                extras[f"pto_{pto.name}_velocity"],
                extras[f"pto_{pto.name}_force"],
                extras[f"pto_{pto.name}_absorbed_power"],
                (DirectDriveHistory(*(
                    extras[f"pto_{pto.name}_drive_{field}"] for field in (
                        "shaft_velocity", "shaft_torque", "inertia_torque",
                        "friction_torque", "generator_torque", "current",
                        "voltage", "resistance_loss", "electrical_power",
                        "mechanical_power",
                    )
                )) if isinstance(pto, LinearPTO) and pto.direct_drive is not None
                 else None),
                (LinearGeneratorHistory(*(
                    extras[f"pto_{pto.name}_generator_{key}"] for key in (
                        "flux_d", "flux_q", "angle", "friction_force",
                        "electrical_power", "phase_current", "phase_voltage",
                    )
                )) if isinstance(pto, LinearPTO) and pto.linear_generator is not None
                 else None),
            )
            for pto in (*self.ptos, *self.rotational_ptos)
        }
        flexible_modes = ({
            self._floating_gbm_body.name: FlexibleModeHistory(
                extras["flex_position"], extras["flex_velocity"],
                extras["flex_acceleration"],
            )
        } if self._floating_gbm_body is not None else {})
        if self._floating_gbm_orifice is not None:
            name = self._floating_gbm_orifice["name"]
            ptos[name] = PTOHistory(
                extras["flex_position"][:, 0],
                extras["flex_velocity"][:, 0],
                extras["orifice_force"], extras["orifice_power"],
            )
            coordinates["pitch_unwrapped"] = MotionHistory(
                extras["pitch_unwrapped"],
                response.body_velocity[:, 0, 4],
            )
        body_forces = {
            body.name: extras[f"morison_force_{body.name}"]
            for body in self.bodies
            if body.fixed or any(attached is body
                                 for attached, _ in self._morison_elements)
        }
        return WECResult(response.time, bodies, coordinates, ptos,
                         response.wave_elevation, case, response,
                         flexible_modes, body_forces)

    @staticmethod
    def _pto_case(pto: LinearPTO) -> dict:
        def endpoint(point):
            if isinstance(point, BodyPoint):
                return {"body": point.body.name, "point": point.coordinates()}
            return {"ground": point.coordinates()}

        item = {
            "name": pto.name,
            "kind": "linear_actuator",
            "from": endpoint(pto.from_point),
            "to": endpoint(pto.to_point),
            "damping": (pto.control.gain if pto.control is not None
                        else pto.damping),
            "stiffness": pto.stiffness,
        }
        if pto.control is not None:
            if isinstance(pto.control, DeclutchingControl):
                item["control"] = {
                    "kind": "declutching",
                    "declutch_time": pto.control.declutch_time,
                    "minimum_on_time": pto.control.minimum_on_time,
                }
            else:
                item["control"] = {
                    "kind": "latching",
                    "latch_time": pto.control.latch_time,
                    "latch_damping": pto.control.latch_damping,
                    "minimum_normal_time": pto.control.minimum_normal_time,
                }
        if pto.axis is not None:
            item["axis"] = list(pto.axis)
        if pto.equilibrium_position is not None:
            item["equilibrium_position"] = pto.equilibrium_position
        if pto.pretension is not None:
            item["pretension"] = pto.pretension
        if pto.direct_drive is not None:
            item["direct_drive"] = vars(pto.direct_drive).copy()
        if pto.linear_generator is not None:
            item["linear_generator"] = vars(pto.linear_generator).copy()
        return item

    @staticmethod
    def _rotational_pto_case(pto: RotationalPTO) -> dict:
        return {
            "name": pto.name,
            "kind": "coordinate_torque",
            "coordinate": pto.coordinate.name,
            "damping": pto.damping,
            "stiffness": pto.stiffness,
            "equilibrium_position": pto.equilibrium_angle,
        }


def _state(value: Mapping[str, float] | Sequence[float]) -> dict | list:
    if isinstance(value, Mapping):
        return dict(value)
    return list(value)

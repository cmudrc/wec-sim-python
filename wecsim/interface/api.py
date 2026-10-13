"""Python builder for WEC configuration and simulation."""

from __future__ import annotations

from numbers import Real
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np

from wecsim.components.directLinearGenerator import DirectLinearGenerator
from wecsim.components.floatingOwc import FloatingOwcChamber, FloatingOwcTurbine
from wecsim.components.orifice import OrificePTO
from wecsim.interface.configuration import (
    Body,
    BodyPoint,
    Coordinate,
    FullDirectionalSpectrumWave,
    HydroState,
    ImportedElevationWave,
    ImportedSpectrumWave,
    JONSWAPWave,
    LinearPTO,
    Motion,
    MultiPMWave,
    NoWave,
    PMWave,
    RegularCICWave,
    RegularWave,
    RotationalPTO,
    SimpleDirectDrive,
    VariableHydro,
    WorldPoint,
    _FixedHinge,
    _FloatingJoint,
    _FloatingOwc,
)
from wecsim.control.controls import DeclutchingControl, LatchingControl
from wecsim.dynamics.caseDynamics import run_case
from wecsim.dynamics.hardStops import LinearHardStops
from wecsim.dynamics.moorDyn import MoorDyn
from wecsim.hydrodynamics.morison import MorisonElement
from wecsim.interface.results import (
    DirectDriveHistory,
    FlexibleModeHistory,
    LinearGeneratorHistory,
    MotionHistory,
    PTOHistory,
    WECResult,
)
from wecsim.interface.wec_case import _WECCaseMixin


class WEC(_WECCaseMixin):
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
                        phase_mode: str = "directional",
                        axis: Sequence[float] | None = None) -> MorisonElement:
        """Attach a Morison element at a body-local point.

        A body-local ``axis`` selects normal/tangential option 2 for a moving
        regular-wave surge/heave/pitch body. Its first two coefficient and
        area entries are normal and tangential. Without an axis, moving
        hydrodynamic bodies support axial Cartesian elements.
        """
        if not any(body is item for item in self.bodies):
            raise ValueError("Morison body must belong to this WEC")
        if not isinstance(point, BodyPoint) or point.body is not body:
            raise ValueError("Morison point must belong to its body")
        if axis is not None:
            vector = np.asarray(axis, dtype=float)
            if (vector.shape != (3,) or not np.isfinite(vector).all()
                    or np.linalg.norm(vector) == 0):
                raise ValueError("Morison axis needs a finite nonzero body-local vector")
        element = MorisonElement(
            tuple(point.coordinates()), tuple(drag_coefficient),
            tuple(added_mass_coefficient), tuple(area), volume, phase_mode,
            tuple(float(value) for value in axis) if axis is not None else None,
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
                       equilibrium_position: float | None = None,
                       pretension: float | None = None,
                       mooring_surge_stiffness: float = 0.0,
                       moordyn: MoorDyn | None = None,
                       moordyn_point: BodyPoint | None = None,
                       hard_stops: LinearHardStops | None = None,
                       radiation_method: str | None = None,
                       added_mass_scheme: str = "implicit") -> None:
        """Select the paired two-body surge/heave/pitch slider joint.

        The PTO acts on float heave minus spar heave. The reduced joint does
        not model off-axis PTO endpoints or arbitrary Simscape constraints.
        Pretension is in newtons and sets equilibrium to -pretension/stiffness.
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
        if radiation_method not in (None, "constant", "convolution", "fir",
                                    "state_space"):
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
        equilibrium = _pto_equilibrium(stiffness, equilibrium_position,
                                       pretension, "floating_joint")
        self._floating_joint = _FloatingJoint(
            float_body, spar_body, location, pto_name, damping, stiffness,
            equilibrium, mooring_surge_stiffness, hard_stops,
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
                    stiffness: float = 0.0,
                    equilibrium_angle: float | None = None,
                    pretension: float | None = None,
                    added_mass_scheme: str = "implicit") -> None:
        """Select the published OSWEC pitch hinge and torsional PTO.

        The flap moves in surge, heave, and pitch about a fixed base. Both
        joint and PTO points currently lie on the world z axis. Pretension
        is in newton meters and sets equilibrium to -pretension/stiffness
        radians when stiffness is in newton meters per radian.
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
        equilibrium = _pto_equilibrium(stiffness, equilibrium_angle,
                                       pretension, "fixed_hinge")
        if (not isinstance(damping, Real) or isinstance(damping, bool)
                or not np.isfinite(damping) or damping < 0):
            raise ValueError("fixed_hinge PTO needs finite nonnegative damping and stiffness")
        if added_mass_scheme not in ("implicit", "simulink_delay"):
            raise ValueError("unsupported fixed_hinge added-mass scheme")
        self._fixed_hinge = _FixedHinge(
            flap, base, location, pto_location, pto_name,
            damping, stiffness, equilibrium, added_mass_scheme,
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
            linear_generator: DirectLinearGenerator | None = None,
            velocity_control_gain: float = 0.0) -> LinearPTO:
        if any(existing.name == name for existing in
               (*self.ptos, *self.rotational_ptos)):
            raise ValueError(f"PTO name already exists: {name}")
        for point in (from_point, to_point):
            if not isinstance(point, (BodyPoint, WorldPoint)):
                raise TypeError("PTO endpoints must be body points or world points")
            if (isinstance(point, BodyPoint)
                    and not any(point.body is body for body in self.bodies)):
                raise ValueError("PTO body point must belong to this WEC")
        try:
            gain = (float(velocity_control_gain)
                    if not isinstance(velocity_control_gain, bool) else np.nan)
        except (TypeError, ValueError) as exc:
            raise ValueError("velocity_control_gain must be finite and nonnegative") from exc
        if not np.isfinite(gain) or gain < 0:
            raise ValueError("velocity_control_gain must be finite and nonnegative")
        if control is not None:
            if not isinstance(control, (DeclutchingControl, LatchingControl)):
                raise TypeError("PTO control must be a supported control law")
            if (damping or stiffness or equilibrium_position is not None
                    or pretension is not None or gain):
                raise ValueError("sampled PTO control sets its own gain and has no spring")
        if direct_drive is not None:
            if not isinstance(direct_drive, SimpleDirectDrive):
                raise TypeError("direct_drive must be SimpleDirectDrive")
            if (control is not None or damping or stiffness or gain
                    or equilibrium_position is not None or pretension is not None
                    or linear_generator is not None):
                raise ValueError("direct drive supplies its own PTO force")
        if linear_generator is not None:
            if not isinstance(linear_generator, DirectLinearGenerator):
                raise TypeError("linear_generator must be DirectLinearGenerator")
            if (control is not None or damping or stiffness or gain
                    or equilibrium_position is not None or pretension is not None):
                raise ValueError("linear generator supplies its own PTO force")
        pto = LinearPTO(
            name, from_point, to_point,
            tuple(axis) if axis is not None else None,
            damping, stiffness, equilibrium_position, pretension, control,
            direct_drive, linear_generator, gain,
        )
        self.ptos.append(pto)
        return pto


    def rotational_pto(self, name: str, coordinate: Coordinate, *,
                       damping: float = 0.0, stiffness: float = 0.0,
                       equilibrium_angle: float | None = None,
                       pretension: float | None = None) -> RotationalPTO:
        """Attach a torsional PTO to a named rotation coordinate.

        Damping is in N m s/rad, stiffness in N m/rad, and equilibrium angle
        in radians. Physical pretension is in N m and sets equilibrium to
        ``-pretension / stiffness``. The returned history reports angle as
        ``stroke``.
        """
        if not any(coordinate is item for item in self.coordinates):
            raise ValueError("rotational PTO coordinate must belong to this WEC")
        if any(existing.name == name for existing in
               (*self.ptos, *self.rotational_ptos)):
            raise ValueError(f"PTO name already exists: {name}")
        equilibrium = _pto_equilibrium(stiffness, equilibrium_angle,
                                       pretension, "rotational")
        pto = RotationalPTO(name, coordinate, damping, stiffness, equilibrium)
        self.rotational_ptos.append(pto)
        return pto


    def run(
        self, wave: RegularWave | RegularCICWave | PMWave | MultiPMWave | JONSWAPWave | ImportedSpectrumWave | FullDirectionalSpectrumWave | ImportedElevationWave | NoWave, *,
        dt: float, end_time: float,
        ramp_time: float | None = None,
        radiation_memory: float | None = None,
        solver: str = "implicit",
        rho: float | None = None, g: float | None = None,
        initial_coordinate: Mapping[str, float] | Sequence[float] | None = None,
        initial_speed: Mapping[str, float] | Sequence[float] | None = None,
        base_dir: str | Path = ".",
    ) -> WECResult:
        """Simulate the WEC and return named NumPy motion and PTO histories.

        ``wecsim_ode4`` is an opt-in source solver for PM sampled yaw;
        the default remains the implicit physical radiation step.
        """
        case = self.to_case(
            wave, dt=dt, end_time=end_time, ramp_time=ramp_time,
            radiation_memory=radiation_memory, solver=solver, rho=rho, g=g,
            initial_coordinate=initial_coordinate, initial_speed=initial_speed,
        )
        # The Python builder and file-driven CLI share the validated case
        # runner; this layer only restores names and structured histories.
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
            # The hinge PTO is rotational, so its stroke and speed are the
            # flap's pitch angle and angular velocity.
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
                (extras[f"pto_{pto.name}_control_force"]
                 if isinstance(pto, LinearPTO) and pto.velocity_control_gain
                 else None),
                (extras[f"pto_{pto.name}_control_absorbed_power"]
                 if isinstance(pto, LinearPTO) and pto.velocity_control_gain
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


def _pto_equilibrium(stiffness: float, equilibrium: float | None,
                     pretension: float | None, layout: str) -> float:
    if equilibrium is not None and pretension is not None:
        raise ValueError(f"{layout} PTO accepts either equilibrium or pretension")
    for name, value in (("stiffness", stiffness), ("equilibrium", equilibrium),
                        ("pretension", pretension)):
        if value is not None and (not isinstance(value, Real)
                                  or isinstance(value, bool)
                                  or not np.isfinite(value)):
            raise ValueError(f"{layout} PTO {name} must be finite")
    if stiffness < 0:
        raise ValueError(f"{layout} PTO stiffness must be nonnegative")
    if (equilibrium or pretension) and not stiffness:
        raise ValueError(
            f"{layout} PTO nonzero equilibrium or pretension needs positive stiffness"
        )
    # Force = -k*(stroke - equilibrium). A positive requested pretension
    # therefore corresponds to a negative equilibrium stroke.
    offset = (-pretension / stiffness if pretension is not None and stiffness
              else equilibrium if equilibrium is not None else 0.0)
    if not np.isfinite(offset):
        raise ValueError(f"{layout} PTO equilibrium is outside its supported range")
    return float(offset)

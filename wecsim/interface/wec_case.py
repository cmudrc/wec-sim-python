"""WEC builder to case serialization."""

from __future__ import annotations

from dataclasses import asdict
from numbers import Real
from typing import Mapping, Sequence

import numpy as np

from wecsim.interface.configuration import (
    BodyPoint,
    FullDirectionalSpectrumWave,
    ImportedElevationWave,
    ImportedSpectrumWave,
    JONSWAPWave,
    LinearPTO,
    MultiPMWave,
    NoWave,
    PMWave,
    RegularCICWave,
    RegularWave,
    RotationalPTO,
    WorldPoint,
)
from wecsim.control.controls import DeclutchingControl


class _WECCaseMixin:
    """Translate a Python WEC builder into a validated case."""

    def to_case(
        self, wave: RegularWave | RegularCICWave | PMWave | MultiPMWave | JONSWAPWave | ImportedSpectrumWave | FullDirectionalSpectrumWave | ImportedElevationWave | NoWave, *,
        dt: float, end_time: float,
        ramp_time: float | None = None,
        radiation_memory: float | None = None,
        solver: str = "implicit",
        rho: float | None = None, g: float | None = None,
        initial_coordinate: Mapping[str, float] | Sequence[float] | None = None,
        initial_speed: Mapping[str, float] | Sequence[float] | None = None,
    ) -> dict:
        """Return the case mapping used by the validated dynamics runner."""
        if not isinstance(wave, (RegularWave, RegularCICWave, PMWave, MultiPMWave,
                                 JONSWAPWave, ImportedSpectrumWave,
                                 FullDirectionalSpectrumWave,
                                 ImportedElevationWave, NoWave)):
            raise TypeError("wave must be a supported regular, irregular, imported, or no-wave configuration")
        if isinstance(wave, MultiPMWave) and self._fixed_hinge is None:
            raise ValueError("MultiPMWave currently needs a fixed-hinge WEC")
        simulation = {"dt": dt, "end_time": end_time}
        if solver != "implicit":
            simulation["solver"] = solver
        for key, value in (
            ("ramp_time", ramp_time), ("radiation_memory", radiation_memory),
            ("rho", rho), ("g", g),
        ):
            if value is not None:
                simulation[key] = value
        # These joints have their own kinematics; only the remaining layouts
        # use the generic linear coordinate maps assembled below.
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
        # Keep builder order: an omitted hydro_body means the matching
        # one-based body index in the shared HDF5 data.
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
                         "phase_mode": element.phase_mode,
                         **({"axis": list(element.axis)}
                            if element.axis is not None else {})}
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
                     "phase_mode": element.phase_mode,
                     **({"axis": list(element.axis)}
                        if element.axis is not None else {})}
                    for element in elements
                ]
            bodies.append(body_case)
        constraint = {"kind": "linear_subspace", "coordinates": []}
        # Fixed nonhydrodynamic bodies have no generalized coordinates;
        # their Morison forces are evaluated on the prescribed pose.
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
        # A named coordinate may contribute motion to several bodies. The
        # case adapter turns these entries into each body's motion Jacobian.
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
        if not isinstance(wave, (RegularWave, PMWave, MultiPMWave,
                                 FullDirectionalSpectrumWave)):
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


    @staticmethod
    def _pto_case(pto: LinearPTO) -> dict:
        def endpoint(point):
            # Retain the body identity and local offset here. The solver uses
            # the HDF5 center of gravity to form the world attachment arm.
            if isinstance(point, BodyPoint):
                return {"body": point.body.name, "point": point.coordinates()}
            return {"ground": point.coordinates()}

        item = {
            "name": pto.name,
            "kind": "linear_actuator",
            "from": endpoint(pto.from_point),
            "to": endpoint(pto.to_point),
            # A sampled controller supplies its own damping gain; the solver
            # advances its switching memory separately from passive damping.
            "damping": (pto.control.gain if pto.control is not None
                        else pto.damping),
            "stiffness": pto.stiffness,
        }
        if pto.velocity_control_gain:
            item["velocity_control_gain"] = pto.velocity_control_gain
        if pto.position_control_gain:
            item["position_control_gain"] = pto.position_control_gain
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

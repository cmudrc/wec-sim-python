"""Specialized heave and floating OWC case adapters."""

import numpy as np

from wecsim.components.directDrive import integrate_direct_drive_heave
from wecsim.components.floatingOwc import FloatingOwcChamber, FloatingOwcTurbine
from wecsim.dynamics.case_inputs import _body_number, _number, _section
from wecsim.dynamics.case_response import CaseResponse
from wecsim.dynamics.floatingOwcDynamics import solve_floating_owc
from wecsim.dynamics.linearCoordinates import build_coordinate_maps, initial_coordinate
from wecsim.dynamics.moorDyn import MoorDyn
from wecsim.dynamics.ptoConnections import build_linear_ptos
from wecsim.hydrodynamics.bodyClass import BodyClass
from wecsim.hydrodynamics.variableHydro import HeaveHydroState, integrate_variable_heave


def _run_floating_owc(case, sim, wave, constraint, bodies, hydro,
                      b2b, dt, end_time, ramp_time, rho, g):
    """Adapt the coupled OWC solver to named case and Python-builder output."""
    if (len(bodies) != 2 or b2b or hydro[0] != hydro[1]
            or "ptos" in case or "radiation_memory" in sim
            or "radiation_method" in sim or "added_mass_scheme" in sim):
        raise ValueError("floating_owc needs one ordered two-body HDF5 and fixed-frequency radiation")
    for index, body in enumerate(bodies, start=1):
        _section(body, "floating_owc body",
                 {"hydro_file", "hydro_body", "mass", "inertia"},
                 {"name", "hydro_file", "hydro_body", "mass", "inertia"})
        _body_number(body, index)
    if (wave["type"] != "regular"
            or set(wave) - {"type", "height", "period", "direction"}
            or _number(wave.get("direction", 0), "wave.direction") != 0):
        raise ValueError("floating_owc needs zero-heading regular waves")
    height = _number(wave.get("height"), "wave.height", nonnegative=True)
    period = _number(wave.get("period"), "wave.period", positive=True)
    _section(constraint, "floating_owc constraint",
             {"kind", "column_height", "column_diameter"},
             {"kind", "column_height", "column_diameter",
              "initial_rotor_speed", "chamber", "turbine"})
    chamber = constraint.get("chamber")
    turbine = constraint.get("turbine")
    if (chamber is not None and not isinstance(chamber, FloatingOwcChamber)
            or turbine is not None and not isinstance(turbine, FloatingOwcTurbine)):
        raise TypeError("floating_owc chamber and turbine need their public classes")
    pto = _section(case.get("pto"), "floating_owc PTO",
                   {"kind", "name", "stiffness", "damping"},
                   {"kind", "name", "stiffness", "damping"})
    if pto["kind"] != "column_slider" or not isinstance(pto["name"], str) or not pto["name"]:
        raise ValueError("floating_owc needs a named column_slider PTO")
    mooring = _section(case.get("mooring"), "floating_owc mooring",
                       {"kind", "session"}, {"kind", "session", "point"})
    if mooring["kind"] != "moor_dyn" or not isinstance(mooring["session"], MoorDyn):
        raise ValueError("floating_owc needs a MoorDyn session")
    point = mooring.get("point")
    if point is not None:
        point = np.asarray(point, dtype=float)
        if point.shape != (3,) or not np.isfinite(point).all():
            raise ValueError("floating_owc mooring point needs three finite local values")
    inertia = [np.asarray(body["inertia"], dtype=float) for body in bodies]
    if (any(value.shape != (3,) or not np.isfinite(value).all()
            for value in inertia) or np.any(inertia[0] <= 0)
            or (np.any(inertia[1] != 0) and np.any(inertia[1] <= 0))):
        raise ValueError("floating_owc body inertias are invalid")
    column_mass = _number(bodies[1]["mass"], "column.mass", positive=True)
    floater_mass = bodies[0]["mass"]
    if floater_mass != "equilibrium":
        floater_mass = _number(floater_mass, "floater.mass", positive=True)
    owc_ramp_time = _number(sim.get("ramp_time", 50.0),
                            "simulation.ramp_time", nonnegative=True)
    response = solve_floating_owc(
        hydro[0], mooring["session"], dt=dt, end_time=end_time,
        wave_height=height, wave_period=period,
        ramp_time=owc_ramp_time,
        floater_mass=floater_mass, floater_inertia=tuple(inertia[0]),
        column_mass=column_mass,
        column_height=_number(constraint["column_height"], "column_height", positive=True),
        column_diameter=_number(constraint["column_diameter"], "column_diameter", positive=True),
        column_inertia=(None if np.all(inertia[1] == 0) else tuple(inertia[1])),
        pto_stiffness=_number(pto["stiffness"], "pto.stiffness"),
        pto_damping=_number(pto["damping"], "pto.damping"),
        moordyn_point=point, chamber=chamber, turbine=turbine,
        initial_rotor_speed=_number(constraint.get("initial_rotor_speed", 150),
                                    "initial_rotor_speed", positive=True),
        rho=rho, g=g,
    )
    wave_ramp = np.ones_like(response.time)
    if owc_ramp_time > 0:
        early = response.time < owc_ramp_time
        wave_ramp[early] = (1 - np.cos(np.pi * response.time[early] / owc_ramp_time)) / 2
    elevation = height / 2 * wave_ramp * np.cos(2 * np.pi * response.time / period)
    return CaseResponse(
        response.time,
        np.stack((response.floater_pose, response.column_pose), axis=1),
        np.stack((response.floater_velocity, response.column_velocity), axis=1),
        hydro, pto_force=response.pto_force,
        pto_label="pto_column_slider_force", wave_elevation=elevation,
        extra_outputs=(
            ("pto_mechanical_power", response.pto_mechanical_power),
            ("chamber_pressure", response.chamber_pressure),
            ("turbine_speed", response.turbine_speed),
            ("turbine_power", response.turbine_power),
            ("pneumatic_power", response.pneumatic_power),
            ("moordyn_connection_position", response.mooring_pose),
            ("moordyn_connection_velocity", response.mooring_velocity),
            ("moordyn_connection_force", response.mooring_force),
        ),
        pto_stroke=response.stroke, pto_velocity=response.stroke_speed,
        pto_absorbed_power=response.pto_dissipated_power,
    )


def _run_direct_drive_heave(case, sim, wave, constraint, bodies, hydro,
                            b2b, dt, end_time, ramp_time, rho, g):
    """Run a regular-wave heave body with the simple direct-drive PTO block."""
    heave_map = np.zeros((6, 1))
    heave_map[2, 0] = 1
    if (len(bodies) != 1 or b2b or wave["type"] != "regular"
            or _number(wave.get("direction", 0), "wave.direction") != 0
            or set(wave) - {"type", "height", "period", "direction"}
            or set(constraint) - {"kind", "coordinates", "initial_coordinate",
                                  "initial_speed"}
            or "radiation_memory" in sim or "pto" in case
            or "mooring" in case or len(case["ptos"]) != 1):
        raise ValueError("simple direct drive needs one heave body, regular zero-heading waves, and one PTO")
    spec = bodies[0]
    if set(spec) - {"name", "hydro_file", "hydro_body", "mass", "inertia"}:
        raise ValueError("direct-drive heave body has unsupported settings")
    _body_number(spec, 1)
    body = BodyClass(str(hydro[0]))
    body.bodyNumber = body.bodyTotal = 1
    body.readH5file()
    if int(np.asarray(body.dof).item()) != 6:
        raise ValueError("direct-drive heave needs a six-DOF BEM dataset")
    center = np.asarray(body.cg, dtype=float).ravel()
    if center.shape != (3,) or not np.isfinite(center).all():
        raise ValueError("direct-drive body needs a finite center")
    name = spec.get("name", "body1")
    maps, coordinate_names = build_coordinate_maps(
        constraint, bodies, [name], [center],
    )
    if len(maps) != 1 or not np.array_equal(maps[0], heave_map):
        raise ValueError("simple direct drive currently supports pure heave")
    connections, _, _, _ = build_linear_ptos(
        case["ptos"], maps, [center], [name], coordinate_names,
        allow_direct_drive=True,
    )
    connection = connections[0]
    if (not np.isclose(connection.stroke_jacobian[0], 1, rtol=0, atol=1e-12)
            or "direct_drive" not in case["ptos"][0]):
        raise ValueError("simple direct drive needs a vertical PTO with positive heave stroke")
    drive = _section(case["ptos"][0]["direct_drive"], "pto.direct_drive",
                     {"kp", "ki", "torque_constant", "gear_ratio",
                      "drivetrain_inertia", "drivetrain_friction",
                      "winding_resistance", "winding_inductance"},
                     {"kp", "ki", "torque_constant", "gear_ratio",
                      "drivetrain_inertia", "drivetrain_friction",
                      "winding_resistance", "winding_inductance"})
    drive = {
        key: _number(value, f"pto.direct_drive.{key}",
                     positive=key in ("torque_constant", "gear_ratio",
                                      "winding_resistance", "winding_inductance"),
                     nonnegative=key in ("drivetrain_inertia", "drivetrain_friction"))
        for key, value in drive.items()
    }
    height = _number(wave.get("height"), "wave.height", nonnegative=True)
    period = _number(wave.get("period"), "wave.period", positive=True)
    mass = spec.get("mass", "equilibrium")
    mass = (rho * float(np.asarray(body.dispVol).item()) if mass == "equilibrium"
            else _number(mass, "body.mass", positive=True))
    inertia = np.asarray(spec.get("inertia", [0, 0, 0]), dtype=float)
    if (inertia.shape != (3,) or not np.isfinite(inertia).all()
            or np.any(inertia < 0)):
        raise ValueError("body.inertia needs three nonnegative values")
    body.mass = mass
    body.hydroStiffness = np.zeros((6, 6))
    body.viscDrag = {"Drag": np.zeros((6, 6)), "cd": np.zeros(6),
                     "characteristicArea": np.zeros(6)}
    body.linearDamping = np.zeros((6, 6))
    frequency = 2 * np.pi / period
    body.hydroForcePre(
        frequency, [0], 1, np.array([0.0]), [], dt, rho, g,
        "regular", np.zeros((2, 1)), 1, 1, 0, 0, 0,
    )
    force = body.hydroForce
    initial_q = initial_coordinate(
        constraint.get("initial_coordinate", [0]), coordinate_names,
        "constraint.initial_coordinate",
    )
    initial_v = initial_coordinate(
        constraint.get("initial_speed", [0]), coordinate_names,
        "constraint.initial_speed",
    )
    response = integrate_direct_drive_heave(
        dt=dt, end_time=end_time, height=height, period=period,
        ramp_time=ramp_time, rho=rho, g=g, mass=mass,
        displaced_volume=float(np.asarray(body.dispVol).item()),
        center_z=float(center[2]),
        added_mass=float(force["fAddedMass"][2, 2]),
        radiation_damping=float(force["fDamping"][2, 2]),
        hydrostatic_stiffness=float(force["linearHydroRestCoef"][2, 2]),
        excitation_real=float(force["fExt"]["re"][2]),
        excitation_imaginary=float(force["fExt"]["im"][2]),
        initial_position=float(center[2]) + initial_q[0],
        initial_velocity=initial_v[0], **drive,
    )
    coordinates = (response.position - center[2])[:, None]
    speeds = response.velocity[:, None]
    positions = np.zeros((len(response.time), 1, 6))
    velocities = np.zeros_like(positions)
    positions[:, 0, 2] = response.position
    velocities[:, 0, 2] = response.velocity
    ramp = np.ones(len(response.time))
    if ramp_time:
        early = response.time < ramp_time
        ramp[early] = (1 - np.cos(np.pi * response.time[early] / ramp_time)) / 2
    elevation = height / 2 * ramp * np.cos(frequency * response.time)
    drive_fields = (
        "shaft_velocity", "shaft_torque", "inertia_torque",
        "friction_torque", "generator_torque", "current", "voltage",
        "resistance_loss", "electrical_power", "mechanical_power",
    )
    extras = (
        (f"coordinate_{coordinate_names[0]}_position", coordinates[:, 0]),
        (f"coordinate_{coordinate_names[0]}_velocity", speeds[:, 0]),
        ("body1_controller_force", response.controller_force),
    ) + tuple(
        (f"pto_{connection.name}_drive_{field}", getattr(response, field))
        for field in drive_fields
    ) + connection.outputs(
        coordinates, speeds, force_override=response.body_force,
    )
    return CaseResponse(
        response.time, positions, velocities, hydro,
        wave_elevation=elevation,
        pto_generalized_force=response.body_force[:, None],
        extra_outputs=extras,
    )


def _run_variable_heave(case, sim, wave, constraint, bodies, hydro,
                        b2b, dt, end_time, ramp_time, rho, g, base_dir):
    """Run the published variable-draft mode with one active heave coordinate."""
    heave_map = np.zeros((6, 1))
    heave_map[2, 0] = 1
    if (len(bodies) != 1 or b2b or wave["type"] != "regular"
            or _number(wave.get("direction", 0), "wave.direction") != 0
            or set(wave) - {"type", "height", "period", "direction"}
            or set(constraint) - {"kind", "coordinates", "initial_coordinate",
                                  "initial_speed"}
            or "radiation_memory" in sim or "pto" in case
            or "mooring" in case):
        raise ValueError("variable hydro currently needs one heave body, regular zero-heading waves, and linear PTO connections")
    spec = bodies[0]
    if (set(spec) - {"name", "hydro_file", "hydro_body", "mass",
                     "inertia", "variable_hydro"}
            or "variable_hydro" not in spec):
        raise ValueError("variable heave body only uses its ordered hydro states")
    _body_number(spec, 1)
    variable = _section(spec["variable_hydro"], "body.variable_hydro",
                        {"states", "switch_times"}, {"states", "switch_times"})
    raw_states = variable["states"]
    if not isinstance(raw_states, list) or len(raw_states) < 2:
        raise ValueError("body.variable_hydro.states needs at least two datasets")
    if not isinstance(variable["switch_times"], list):
        raise ValueError("body.variable_hydro.switch_times must be a list")
    paths = []
    for number, entry in enumerate(raw_states, start=1):
        entry = _section(entry, f"body.variable_hydro.states[{number}]",
                         {"hydro_file", "mass"},
                         {"hydro_file", "mass", "inertia"})
        raw = entry["hydro_file"]
        if not isinstance(raw, str) or not raw:
            raise ValueError("hydro state file must be a path")
        path = (base_dir / raw).expanduser().resolve(strict=True)
        if not path.is_file():
            raise ValueError(f"hydro state file is not a file: {path}")
        inertia = np.asarray(entry.get("inertia", [0, 0, 0]), dtype=float)
        if (inertia.shape != (3,) or not np.isfinite(inertia).all()
                or np.any(inertia < 0)):
            raise ValueError("hydro state inertia needs three nonnegative values")
        paths.append(path)
    if paths[0] != hydro[0]:
        raise ValueError("body.hydro_file must be the first variable hydro state")

    height = _number(wave.get("height"), "wave.height", nonnegative=True)
    period = _number(wave.get("period"), "wave.period", positive=True)
    frequency = 2 * np.pi / period
    loaded = []
    for path in paths:
        body = BodyClass(str(path))
        body.bodyNumber = 1
        body.bodyTotal = 1
        body.readH5file()
        if int(np.asarray(body.dof).item()) != 6:
            raise ValueError("variable heave needs six-DOF BEM datasets")
        loaded.append(body)
    centers = [np.asarray(loaded[0].cg, dtype=float).ravel()]
    if centers[0].shape != (3,) or not np.isfinite(centers[0]).all():
        raise ValueError("variable heave needs a finite initial center")
    names = [spec.get("name", "body1")]
    maps, coordinate_names = build_coordinate_maps(
        constraint, bodies, names, centers,
    )
    if len(maps) != 1 or not np.array_equal(maps[0], heave_map):
        raise ValueError("variable hydro currently supports pure heave motion")
    if "ptos" not in case or len(case["ptos"]) != 1:
        raise ValueError("variable heave needs one linear PTO connection")
    connections, stiffness, damping, bias = build_linear_ptos(
        case["ptos"], maps, centers, names, coordinate_names,
    )
    connection = connections[0]
    if connection.control is not None or not np.isclose(abs(connection.stroke_jacobian[0]), 1):
        raise ValueError("variable heave needs a passive vertical PTO")
    initial_q = initial_coordinate(
        constraint.get("initial_coordinate", [0]), coordinate_names,
        "constraint.initial_coordinate",
    )
    initial_v = initial_coordinate(
        constraint.get("initial_speed", [0]), coordinate_names,
        "constraint.initial_speed",
    )

    states = []
    for entry, body in zip(raw_states, loaded):
        mass = (rho * float(np.asarray(body.dispVol).item())
                if entry["mass"] == "equilibrium" else
                _number(entry["mass"], "hydro state mass", positive=True))
        body.mass = mass
        body.hydroStiffness = np.zeros((6, 6))
        body.viscDrag = {"Drag": np.zeros((6, 6)), "cd": np.zeros(6),
                         "characteristicArea": np.zeros(6)}
        body.linearDamping = np.zeros((6, 6))
        body.hydroForcePre(
            frequency, [0], 1, np.array([0.0]), [], dt, rho, g,
            "regular", np.zeros((2, 1)), 1, 1, 0, 0, 0,
        )
        force = body.hydroForce
        center = np.asarray(body.cg, dtype=float).ravel()
        if center.shape != (3,) or not np.isfinite(center).all():
            raise ValueError("hydro state center must be finite")
        states.append(HeaveHydroState(
            mass=mass,
            displaced_volume=float(np.asarray(body.dispVol).item()),
            equilibrium_z=float(center[2]),
            added_mass=float(force["fAddedMass"][2, 2]),
            radiation_damping=float(force["fDamping"][2, 2]),
            hydrostatic_stiffness=float(force["linearHydroRestCoef"][2, 2]),
            excitation_real=float(force["fExt"]["re"][2]),
            excitation_imaginary=float(force["fExt"]["im"][2]),
        ))
    response = integrate_variable_heave(
        tuple(states), tuple(variable["switch_times"]),
        dt=dt, end_time=end_time, height=height, period=period,
        ramp_time=ramp_time, rho=rho, g=g,
        pto_stiffness=float(stiffness[0, 0]),
        pto_damping=float(damping[0, 0]), pto_bias=float(bias[0]),
        initial_position=states[0].equilibrium_z + initial_q[0],
        initial_velocity=initial_v[0],
    )
    coordinates = (response.position - states[0].equilibrium_z)[:, None]
    speeds = response.velocity[:, None]
    positions = np.zeros((len(response.time), 1, 6))
    velocities = np.zeros_like(positions)
    positions[:, 0, 2] = response.position
    velocities[:, 0, 2] = response.velocity
    ramp = np.ones(len(response.time))
    if ramp_time:
        early = response.time < ramp_time
        ramp[early] = (1 - np.cos(np.pi * response.time[early] / ramp_time)) / 2
    elevation = height / 2 * ramp * np.cos(frequency * response.time)
    extras = (
        (f"coordinate_{coordinate_names[0]}_position", coordinates[:, 0]),
        (f"coordinate_{coordinate_names[0]}_velocity", speeds[:, 0]),
        ("body1_hydro_state", response.active_state + 1),
        ("body1_excitation_heave", response.excitation),
        ("body1_radiation_heave", response.radiation),
        ("body1_restoring_heave", response.restoring),
    ) + connection.outputs(coordinates, speeds)
    return CaseResponse(
        response.time, positions, velocities, (paths[0],),
        wave_elevation=elevation, auxiliary_files=tuple(paths[1:]),
        pto_generalized_force=response.pto_force[:, None],
        extra_outputs=extras,
    )

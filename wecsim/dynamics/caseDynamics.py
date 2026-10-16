"""Case-driven entry point for supported device layouts.

Parsing and specialized solvers live beside this entry point.
"""

from pathlib import Path
from typing import Mapping

import numpy as np

from wecsim.components.orifice import OrificePTO
from wecsim.dynamics.case_inputs import (
    _body_number,
    _hard_stops,
    _hydro_file,
    _imported_elevation,
    _irregular_components_from_case,
    _location,
    _number,
    _pto,
    _section,
)
from wecsim.dynamics.case_response import CaseResponse
from wecsim.dynamics.gbmFloating import (
    solve_floating_gbm_pm_orifice,
    solve_floating_gbm_regular,
)
from wecsim.dynamics.heave_cases import (
    _run_direct_drive_heave,
    _run_floating_owc,
    _run_variable_heave,
)
from wecsim.dynamics.hingePitch import (
    solve_hinged_pitch_from_excitation,
    solve_hinged_pitch_nonlinear_regular,
    solve_hinged_pitch_regular,
)
from wecsim.dynamics.linearCoordinates import initial_coordinate
from wecsim.dynamics.linearHeave import solve_heave_free_decay
from wecsim.dynamics.morison_case import _run_fixed_morison
from wecsim.dynamics.rm3Regular import solve_rm3_regular
from wecsim.dynamics.subspace_case import _run_linear_subspace
from wecsim.environment.irregularWave import (
    imported_full_directional_components,
    imported_spectrum_components,
    pm_equal_energy_components,
    synthesize_full_directional_response,
    synthesize_irregular_response,
    synthesize_multiple_irregular_response,
)
from wecsim.hydrodynamics.bodyClass import BodyClass

AXES = ("surge", "sway", "heave", "roll", "pitch", "yaw")


def run_case(case: Mapping, *, base_dir: str | Path = ".") -> CaseResponse:
    """Run one explicitly supported wave/device configuration.

    ``base_dir`` resolves relative hydro and phase file paths. Supported
    layouts include heave free decay, a fixed hinge, a floating joint, a
    floating OWC, mapped linear coordinates, and fixed Morison bodies.
    """
    case = _section(
        case, "case", {"simulation", "wave", "bodies", "constraint"},
        {"name", "simulation", "wave", "bodies", "constraint", "pto",
         "ptos", "body_to_body", "mooring"},
    )
    sim = _section(case["simulation"], "simulation", {"dt", "end_time"},
                   {"dt", "end_time", "ramp_time", "rho", "g",
                    "radiation_memory", "radiation_method", "added_mass_scheme",
                    "solver"})
    dt = _number(sim["dt"], "simulation.dt", positive=True)
    end_time = _number(sim["end_time"], "simulation.end_time", nonnegative=True)
    ramp_time = _number(sim.get("ramp_time", 100), "simulation.ramp_time", nonnegative=True)
    rho = _number(sim.get("rho", 1000), "simulation.rho", positive=True)
    g = _number(sim.get("g", 9.81), "simulation.g", positive=True)
    wave = _section(case["wave"], "wave", {"type"},
                    {"type", "height", "period", "direction", "directions",
                     "spreading", "seed", "phase_file", "phase_generator", "frequency_count",
                     "gamma", "file", "variable", "reapply_force_ramp", "seas",
                     "excitation_interpolation", "force_quadrature",
                     "frequency_range", "discretization", "water_depth", "current"})
    constraint = _section(case["constraint"], "constraint", {"kind"},
                          {"kind", "location", "initial_displacement",
                           "initial_coordinate", "initial_speed", "coordinates",
                           "orifice", "orifice_force_path", "heave_linear_damping",
                           "mode_linear_damping", "heave_drag_cd",
                           "heave_drag_area", "pitch_drag_cd",
                           "pitch_drag_area", "column_height",
                           "column_diameter", "initial_rotor_speed",
                           "chamber", "turbine"})
    if "current" in wave and constraint["kind"] not in ("fixed_morison", "linear_subspace"):
        raise ValueError("wave.current requires fixed_morison or moving Morison")
    bodies = case["bodies"]
    if not isinstance(bodies, list) or not bodies:
        raise ValueError("bodies must be a nonempty list")
    base = Path(base_dir).expanduser().resolve()
    hydro = tuple(_hydro_file(body, base) for body in bodies)
    b2b = case.get("body_to_body", False)
    if not isinstance(b2b, bool):
        raise ValueError("body_to_body must be a boolean")
    kind = constraint["kind"]
    for index, body in enumerate(bodies):
        if "fixed" in body and not (
            kind == "fixed_morison" or
            (kind == "fixed_hinge" and len(bodies) == 2
             and index == 1 and body["fixed"] is True)
        ):
            raise ValueError("a fixed hydrodynamic body currently requires a two-body fixed_hinge")
    if "radiation_method" in sim and kind != "floating_joint":
        raise ValueError("radiation_method currently applies to floating_joint")
    if "solver" in sim and kind != "linear_subspace" and not (
            sim["solver"] == "wecsim_ode4" and (
                (kind == "fixed_hinge" and wave["type"] == "regular")
                or (kind == "floating_joint" and wave["type"] == "none"))):
        raise ValueError(
            "simulation.solver needs linear_subspace, a nonlinear regular "
            "fixed_hinge, or a no-wave floating_joint"
        )
    if "added_mass_scheme" in sim and kind not in ("floating_joint", "fixed_hinge"):
        raise ValueError("added_mass_scheme currently applies to floating_joint or fixed_hinge")
    if any(path is None for path in hydro) and not (
        (kind == "fixed_hinge" and len(bodies) == 2
         and hydro[0] is not None and hydro[1] is None)
        or (kind == "fixed_morison" and all(path is None for path in hydro))
    ):
        raise ValueError("a fixed nonhydrodynamic body requires a two-body fixed_hinge")
    if "ptos" in case and kind != "linear_subspace":
        raise ValueError("configurable PTO connections require linear_subspace")
    if kind != "linear_subspace" and any("mean_drift" in body for body in bodies):
        raise ValueError("body.mean_drift currently requires linear_subspace")
    if "mooring" in case and kind not in ("floating_joint", "floating_owc"):
        raise ValueError("the configured mooring requires a floating_joint or floating_owc")

    # Layout-specific adapters retain their own joint or fluid mechanics;
    # mapped coordinates are handled by the shared subspace solver below.
    if kind == "floating_owc":
        return _run_floating_owc(case, sim, wave, constraint, bodies, hydro,
                                 b2b, dt, end_time, ramp_time, rho, g)

    if kind == "fixed_morison":
        return _run_fixed_morison(
            case, sim, wave, constraint, bodies, hydro, b2b,
            dt, end_time, ramp_time, rho, g, base,
        )

    if "ptos" in case and any(
            isinstance(pto, Mapping) and "direct_drive" in pto
            for pto in case["ptos"]):
        if kind != "linear_subspace":
            raise ValueError("simple direct drive needs linear_subspace coordinates")
        return _run_direct_drive_heave(
            case, sim, wave, constraint, bodies, hydro, b2b,
            dt, end_time, ramp_time, rho, g,
        )

    if any("variable_hydro" in body for body in bodies):
        if kind != "linear_subspace":
            raise ValueError("variable hydrodynamics need linear_subspace coordinates")
        return _run_variable_heave(
            case, sim, wave, constraint, bodies, hydro, b2b,
            dt, end_time, ramp_time, rho, g, base,
        )

    if kind == "floating_gbm":
        if (len(bodies) != 1 or b2b
                or "pto" in case or "ptos" in case):
            raise ValueError("floating GBM needs one body and no separate PTO")
        if np.any(_location(constraint)):
            raise ValueError("floating GBM needs a joint at the body origin")
        if set(bodies[0]) - {"hydro_file", "hydro_body", "mass", "inertia", "name"}:
            raise ValueError("floating GBM uses body mass and inertia only")
        _body_number(bodies[0], 1)
        inertia = np.asarray(bodies[0].get("inertia"), dtype=float)
        if inertia.shape != (3,) or not np.isfinite(inertia).all() or np.any(inertia <= 0):
            raise ValueError("floating GBM needs three positive body inertia values")
        body_mass = bodies[0].get("mass", "equilibrium")
        if body_mass != "equilibrium":
            body_mass = _number(body_mass, "body.mass", positive=True)
        auxiliary_files = ()
        if "orifice" not in constraint:
            if (wave["type"] != "regular"
                    or set(sim) - {"dt", "end_time", "ramp_time", "rho", "g"}
                    or set(wave) - {"type", "height", "period", "direction"}
                    or set(constraint) - {"kind", "location"}):
                raise ValueError("floating GBM without an orifice needs one regular wave direction")
            if _number(wave.get("direction", 0), "wave.direction") != 0:
                raise ValueError("floating GBM currently needs a zero-degree wave")
            result = solve_floating_gbm_regular(
                hydro[0], dt=dt, end_time=end_time,
                height=_number(wave.get("height"), "wave.height", nonnegative=True),
                period=_number(wave.get("period"), "wave.period", positive=True),
                pitch_inertia=float(inertia[1]), mass=body_mass,
                ramp_time=ramp_time, rho=rho, g=g,
            )
        else:
            if (wave["type"] != "pm"
                    or set(sim) - {"dt", "end_time", "ramp_time", "rho", "g",
                                   "radiation_memory"}
                    or set(wave) - {"type", "height", "period", "directions",
                                    "spreading", "seed", "phase_file", "phase_generator",
                                    "frequency_count"}
                    or set(constraint) - {"kind", "location", "orifice",
                                          "orifice_force_path",
                                          "heave_linear_damping",
                                          "mode_linear_damping", "heave_drag_cd",
                                          "heave_drag_area", "pitch_drag_cd",
                                          "pitch_drag_area"}):
                raise ValueError("floating OWC orifice needs a zero-heading PM sea")
            directions = np.asarray(wave.get("directions", [0]), dtype=float)
            spreading = np.asarray(wave.get("spreading", [1]), dtype=float)
            if (directions.shape != (1,) or spreading.shape != (1,)
                    or not np.isfinite(directions).all()
                    or not np.isfinite(spreading).all()
                    or directions[0] != 0 or spreading[0] != 1):
                raise ValueError("floating OWC currently needs one zero-heading PM component")
            if "seed" in wave and "phase_file" in wave:
                raise ValueError("PM wave uses either seed or phase_file")
            if "phase_file" in wave:
                phase_path = (base / wave["phase_file"]).resolve(strict=True)
                phase = np.loadtxt(phase_path, delimiter=",", ndmin=2)
                seed = None
                auxiliary_files = (phase_path,)
            else:
                phase = None
                seed = wave.get("seed", 7)
                if not isinstance(seed, int) or isinstance(seed, bool):
                    raise ValueError("wave.seed must be an integer")
            components = pm_equal_energy_components(
                hydro[0],
                significant_height=_number(wave.get("height"), "wave.height", positive=True),
                peak_period=_number(wave.get("period"), "wave.period", positive=True),
                directions=directions, spreading=spreading,
                count=wave.get("frequency_count", 500), seed=seed, phase=phase,
                phase_generator=wave.get("phase_generator", "numpy"),
            )
            orifice_spec = _section(
                constraint["orifice"], "constraint.orifice",
                {"piston_area", "orifice_area"},
                {"name", "piston_area", "orifice_area",
                 "discharge_coefficient", "air_density", "mach_threshold",
                 "sound_speed"},
            )
            orifice = OrificePTO(**{
                key: value for key, value in orifice_spec.items()
                if key != "name"
            })
            result = solve_floating_gbm_pm_orifice(
                hydro[0], dt=dt, end_time=end_time,
                components=components, pitch_inertia=float(inertia[1]),
                orifice=orifice, mass=body_mass, ramp_time=ramp_time,
                orifice_force_path=constraint.get("orifice_force_path", "coupled"),
                memory_time=_number(sim.get("radiation_memory", 15),
                                    "simulation.radiation_memory", positive=True),
                heave_linear_damping=_number(constraint.get("heave_linear_damping", 0),
                                             "constraint.heave_linear_damping", nonnegative=True),
                mode_linear_damping=_number(constraint.get("mode_linear_damping", 0),
                                            "constraint.mode_linear_damping", nonnegative=True),
                heave_drag_cd=_number(constraint.get("heave_drag_cd", 0),
                                      "constraint.heave_drag_cd", nonnegative=True),
                heave_drag_area=_number(constraint.get("heave_drag_area", 0),
                                        "constraint.heave_drag_area", nonnegative=True),
                pitch_drag_cd=_number(constraint.get("pitch_drag_cd", 0),
                                      "constraint.pitch_drag_cd", nonnegative=True),
                pitch_drag_area=_number(constraint.get("pitch_drag_area", 0),
                                        "constraint.pitch_drag_area", nonnegative=True),
                rho=rho, g=g,
            )
        extras = [
            ("flex_position", result.mode_position),
            ("flex_velocity", result.mode_velocity),
            ("flex_acceleration", result.mode_acceleration),
        ]
        if result.orifice_force is not None:
            extras += [
                ("orifice_force", result.orifice_force),
                ("orifice_power", result.orifice_power),
                ("orifice_compressibility_flag",
                 result.orifice_compressibility_flag),
                ("pitch_unwrapped", result.unwrapped_pitch),
            ]
        return CaseResponse(
            result.time, result.body_position[:, None, :],
            result.body_velocity[:, None, :], hydro,
            wave_elevation=result.wave_elevation,
            auxiliary_files=auxiliary_files,
            extra_outputs=tuple(extras),
        )

    if kind == "linear_subspace":
        # Named coordinate maps define the body kinematics and PTO stroke
        # projections without integrating redundant six-DOF body states.
        return _run_linear_subspace(
            case, sim, wave, constraint, bodies, hydro, b2b,
            dt, end_time, ramp_time, rho, g, base,
        )

    if kind == "heave":
        if len(bodies) != 1 or wave["type"] != "none" or b2b or "pto" in case:
            raise ValueError("heave free decay needs one body, no waves, and no PTO")
        if set(wave) != {"type"} or set(constraint) - {"kind", "initial_displacement"}:
            raise ValueError("heave free decay has no wave or joint-location settings")
        if set(bodies[0]) - {"hydro_file", "hydro_body", "mass"}:
            raise ValueError("heave free decay does not use body inertia settings")
        _body_number(bodies[0], 1)
        if bodies[0].get("mass", "equilibrium") != "equilibrium":
            raise ValueError("heave free decay currently requires equilibrium mass")
        if "ramp_time" in sim:
            raise ValueError("ramp_time is inapplicable to no-wave free decay")
        displacement = _number(
            constraint.get("initial_displacement", 0),
            "constraint.initial_displacement",
        )
        solved = solve_heave_free_decay(
            hydro[0], displacement, dt=dt, end_time=end_time,
            cic_end_time=_number(
                sim.get("radiation_memory", 15), "simulation.radiation_memory",
                positive=True,
            ),
            rho=rho, g=g,
        )
        position = np.zeros((len(solved.time), 1, 6))
        velocity = np.zeros_like(position)
        position[:, 0, 2] = solved.position
        velocity[:, 0, 2] = solved.velocity
        return CaseResponse(
            solved.time, position, velocity, hydro,
            total_heave_force=solved.force_total,
        )

    if kind == "fixed_hinge":
        # The hinge solver follows the flap's curved center-of-gravity path;
        # a constant pitch coordinate map would miss that translation.
        if (len(bodies) not in (1, 2)
                or wave["type"] not in ("pm", "jonswap", "pm_multi", "regular",
                                        "spectrumImportFullDir")
                or b2b):
            raise ValueError("fixed-hinge pitch needs one flap, optional fixed base, and regular, PM, JONSWAP, multi-PM, or full-directional waves")
        if len(bodies) == 2 and hydro[1] is not None:
            _body_number(bodies[1], 2)
            fixed_body = BodyClass(str(hydro[1]))
            fixed_body.bodyNumber = 2
            fixed_body.readH5file()
            fixed_center = np.asarray(fixed_body.cg, dtype=float).ravel()
        elif len(bodies) == 2:
            fixed_center = np.asarray(bodies[1]["center_gravity"], dtype=float)
        _body_number(bodies[0], 1)
        if set(bodies[0]) - {"hydro_file", "hydro_body", "mass", "pitch_inertia",
                              "name", "nonlinear_hydro", "geometry_file"}:
            raise ValueError("fixed-hinge pitch uses mass, pitch inertia, and optional mesh hydro")
        nonlinear = bodies[0].get("nonlinear_hydro")
        geometry = bodies[0].get("geometry_file")
        if nonlinear not in (None, "instantaneous") or (geometry is None) != (nonlinear is None):
            raise ValueError("nonlinear fixed hinge needs instantaneous hydro and a mesh")
        if nonlinear is not None and (not isinstance(geometry, str) or not geometry):
            raise ValueError("fixed-hinge geometry_file must be a nonempty file path")
        mass = _number(bodies[0].get("mass"), "body.mass", positive=True)
        inertia = _number(bodies[0].get("pitch_inertia"),
                          "body.pitch_inertia", positive=True)
        if set(constraint) - {"kind", "location", "initial_coordinate", "initial_speed"}:
            raise ValueError("unsupported fixed-hinge constraint setting")
        location = _location(constraint)
        initial_angle = initial_coordinate(
            constraint.get("initial_coordinate", [0]), ("pitch",),
            "constraint.initial_coordinate",
        )[0]
        initial_pitch_speed = initial_coordinate(
            constraint.get("initial_speed", [0]), ("pitch",),
            "constraint.initial_speed",
        )[0]
        pto_data = case.get("pto")
        damping, stiffness, equilibrium = _pto(
            pto_data, "pitch", allow_location=True,
        )
        if len(bodies) == 2 and "location" not in pto_data:
            raise ValueError("fixed nonhydrodynamic base needs pto.location")
        hinge = (_location({"location": pto_data["location"]}, "pto.location")
                 if "location" in pto_data else location)
        if wave["type"] == "regular":
            source_ode4 = sim.get("solver") == "wecsim_ode4"
            scheme = sim.get("added_mass_scheme", "implicit")
            if source_ode4:
                if nonlinear != "instantaneous" or scheme != "simulink_delay":
                    raise ValueError(
                        "regular source ode4 needs nonlinear mesh pressures and simulink_delay"
                    )
            elif scheme != "implicit":
                raise ValueError("regular fixed-hinge dynamics need implicit added mass")
            height = _number(wave.get("height"), "wave.height", nonnegative=True)
            period = _number(wave.get("period"), "wave.period", positive=True)
            if set(wave) - {"type", "height", "period", "direction"}:
                raise ValueError("regular waves use height, period, and direction")
            direction = _number(wave.get("direction", 0), "wave.direction")
            if direction != 0:
                raise ValueError("fixed-hinge regular dynamics currently support 0-degree waves")
            if "radiation_memory" in sim:
                raise ValueError("regular fixed-hinge dynamics use constant-frequency radiation")
            solver = (solve_hinged_pitch_nonlinear_regular if nonlinear is not None
                      else solve_hinged_pitch_regular)
            if nonlinear is not None:
                geometry_path = (base / geometry).resolve(strict=True)
            else:
                geometry_path = None
            source_options = {"source_ode4": True} if source_ode4 else {}
            solved = solver(
                hydro[0], *([geometry_path] if geometry_path is not None else []),
                wave_height=height, wave_period=period,
                hinge_z=hinge[2], body_mass=mass, pitch_inertia=inertia,
                pto_damping=damping, pto_stiffness=stiffness,
                pto_equilibrium=equilibrium, dt=dt, end_time=end_time,
                ramp_time=ramp_time, rho=rho, g=g,
                initial_angle=initial_angle,
                initial_angular_velocity=initial_pitch_speed,
                **source_options,
            )
            ramp = np.ones(len(solved.time))
            if ramp_time > 0:
                early = solved.time < ramp_time
                ramp[early] = (1 - np.cos(np.pi * solved.time[early] / ramp_time)) / 2
            elevation = height / 2 * ramp * np.cos(2 * np.pi * solved.time / period)
            position = np.zeros((len(solved.time), len(bodies), 6))
            velocity = np.zeros_like(position)
            position[:, 0, :3] = solved.center_position
            position[:, 0, 4] = solved.angle
            velocity[:, 0, :3] = solved.center_velocity
            velocity[:, 0, 4] = solved.angular_velocity
            if len(bodies) == 2:
                position[:, 1, :3] = fixed_center
            return CaseResponse(
                solved.time, position, velocity,
                tuple(path for path in hydro if path is not None),
                pto_force=solved.pto_torque, pto_label="pto_pitch_torque",
                wave_elevation=elevation,
            )
        if nonlinear is not None:
            raise ValueError("nonlinear fixed hinge currently needs a regular wave")
        if wave["type"] == "spectrumImportFullDir":
            if (set(wave) - {"type", "file", "phase_file", "seed", "phase_generator",
                             "excitation_interpolation", "force_quadrature"}
                    or "file" not in wave):
                raise ValueError("spectrumImportFullDir needs a MAT spectrum file")
            if "phase_file" in wave and "seed" in wave:
                raise ValueError("supply either wave.phase_file or wave.seed")
            if not isinstance(wave["file"], str) or not wave["file"]:
                raise ValueError("wave.file must be a MAT file path")
            spectrum_path = (base / wave["file"]).resolve(strict=True)
            auxiliary_files = [spectrum_path]
            if "phase_file" in wave:
                if not isinstance(wave["phase_file"], str) or not wave["phase_file"]:
                    raise ValueError("wave.phase_file must be a file path")
                phase_path = (base / wave["phase_file"]).resolve(strict=True)
                phase = np.loadtxt(phase_path, delimiter=",", ndmin=2)
                auxiliary_files.append(phase_path)
                seed = None
            else:
                phase = None
                seed = wave.get("seed", 7)
                if not isinstance(seed, int) or isinstance(seed, bool):
                    raise ValueError("wave.seed must be an integer")
            components = imported_full_directional_components(
                hydro[0], spectrum_path, phase=phase, seed=seed,
                phase_generator=wave.get("phase_generator", "numpy"),
            )
            incident = synthesize_full_directional_response(
                hydro[0], components, dt=dt, end_time=end_time,
                ramp_time=ramp_time, rho=rho, g=g,
                excitation_interpolation=wave.get("excitation_interpolation", "linear"),
                force_quadrature=wave.get("force_quadrature", "integrated"),
            )
        elif wave["type"] == "pm_multi":
            if set(wave) - {"type", "seas", "excitation_interpolation"} or "seas" not in wave:
                raise ValueError("pm_multi waves use a seas list")
            seas = wave["seas"]
            if not isinstance(seas, list) or len(seas) < 2:
                raise ValueError("pm_multi needs at least two sea spectra")
            components = []
            auxiliary_files = []
            for index, sea in enumerate(seas):
                _section(sea, f"wave.seas[{index}]", {"height", "period"},
                         {"height", "period", "direction", "directions",
                          "spreading", "seed", "phase_file", "phase_generator",
                          "frequency_count"})
                if "direction" in sea:
                    if "directions" in sea or "spreading" in sea:
                        raise ValueError("a sea uses direction or directions and spreading")
                    directions = [_number(sea["direction"],
                                          f"wave.seas[{index}].direction")]
                    spreading = [1.0]
                else:
                    if "directions" not in sea or "spreading" not in sea:
                        raise ValueError("a sea needs direction or directions and spreading")
                    directions, spreading = sea["directions"], sea["spreading"]
                if "phase_file" in sea and "seed" in sea:
                    raise ValueError("a sea uses either phase_file or seed")
                if "phase_file" in sea:
                    if not isinstance(sea["phase_file"], str) or not sea["phase_file"]:
                        raise ValueError("sea.phase_file must be a file path")
                    phase_path = (base / sea["phase_file"]).resolve(strict=True)
                    phase = np.loadtxt(phase_path, delimiter=",", ndmin=2)
                    seed = None
                    auxiliary_files.append(phase_path)
                else:
                    phase = None
                    seed = sea.get("seed", index + 7)
                    if not isinstance(seed, int) or isinstance(seed, bool):
                        raise ValueError("sea.seed must be an integer")
                components.append(pm_equal_energy_components(
                    hydro[0],
                    significant_height=_number(
                        sea["height"], f"wave.seas[{index}].height", positive=True,
                    ),
                    peak_period=_number(
                        sea["period"], f"wave.seas[{index}].period", positive=True,
                    ),
                    directions=directions, spreading=spreading,
                    count=sea.get("frequency_count", 500), seed=seed, phase=phase,
                    phase_generator=sea.get("phase_generator", "numpy"),
                ))
            incident = synthesize_multiple_irregular_response(
                hydro[0], components, dt=dt, end_time=end_time,
                ramp_time=ramp_time, rho=rho, g=g,
                excitation_interpolation=wave.get("excitation_interpolation", "linear"),
            )
        else:
            # Preserve the historical fixed-hinge PM spread for case files
            # that omit directions; public PMWave supplies its own settings.
            wave_settings = ({"directions": [0, 30, 90],
                              "spreading": [0.1, 0.2, 0.7], **wave}
                             if wave["type"] == "pm" else wave)
            components, auxiliary_files = _irregular_components_from_case(
                wave_settings, hydro[0], base,
            )
            incident = synthesize_irregular_response(
                hydro[0], components, dt=dt, end_time=end_time,
                ramp_time=ramp_time, rho=rho, g=g,
                excitation_interpolation=wave.get("excitation_interpolation", "linear"),
            )
        solved = solve_hinged_pitch_from_excitation(
            hydro[0], incident.excitation_force,
            hinge_z=hinge[2], body_mass=mass, pitch_inertia=inertia,
            pto_damping=damping, pto_stiffness=stiffness,
            pto_equilibrium=equilibrium, dt=dt,
            memory_time=_number(
                sim.get("radiation_memory", 30), "simulation.radiation_memory",
                positive=True,
            ),
            added_mass_scheme=sim.get("added_mass_scheme", "implicit"),
            rho=rho, g=g,
            initial_angle=initial_angle,
            initial_angular_velocity=initial_pitch_speed,
        )
        position = np.zeros((len(solved.time), len(bodies), 6))
        velocity = np.zeros_like(position)
        position[:, 0, :3] = solved.center_position
        position[:, 0, 4] = solved.angle
        velocity[:, 0, :3] = solved.center_velocity
        velocity[:, 0, 4] = solved.angular_velocity
        if len(bodies) == 2:
            position[:, 1, :3] = fixed_center
        return CaseResponse(
            solved.time, position, velocity,
            tuple(path for path in hydro if path is not None),
            pto_force=solved.pto_torque, pto_label="pto_pitch_torque",
            wave_elevation=incident.elevation,
            auxiliary_files=tuple(auxiliary_files),
        )

    if kind == "floating_joint":
        if len(bodies) != 2 or wave["type"] not in (
                "regular", "regularCIC", "pm", "jonswap", "spectrumImport",
                "none", "elevationImport"):
            raise ValueError("floating joint needs two bodies and supported regular or irregular waves")
        if hydro[0] != hydro[1]:
            raise ValueError("the current floating-joint layout needs one shared HDF5")
        for number, body in enumerate(bodies, start=1):
            _body_number(body, number)
            if set(body) - {"hydro_file", "hydro_body", "mass", "pitch_inertia"}:
                raise ValueError("floating-joint bodies use mass and pitch_inertia")
            if body.get("mass", "equilibrium") != "equilibrium":
                raise ValueError("floating-joint bodies currently require equilibrium mass")
        if set(constraint) - {"kind", "location", "initial_coordinate", "initial_speed"}:
            raise ValueError("unsupported floating-joint constraint setting")
        location = _location(constraint)
        coordinate_names = ("surge", "float_heave", "spar_heave", "pitch")
        initial_q = initial_coordinate(
            constraint.get("initial_coordinate", [0] * 4),
            coordinate_names, "constraint.initial_coordinate",
        )
        initial_v = initial_coordinate(
            constraint.get("initial_speed", [0] * 4),
            coordinate_names, "constraint.initial_speed",
        )
        pto_spec = case.get("pto")
        hard_stops = None
        if isinstance(pto_spec, Mapping) and "hard_stops" in pto_spec:
            hard_stops = _hard_stops(pto_spec["hard_stops"])
            pto_spec = {key: value for key, value in pto_spec.items()
                        if key != "hard_stops"}
        damping, stiffness, equilibrium = _pto(pto_spec, "relative_heave")
        mooring_stiffness = 0.0
        moordyn = moordyn_point = None
        if "mooring" in case:
            mooring = _section(case["mooring"], "mooring", {"kind"},
                               {"kind", "stiffness", "session", "point"})
            if mooring["kind"] == "joint_surge_spring":
                if set(mooring) != {"kind", "stiffness"}:
                    raise ValueError("joint surge spring needs only stiffness")
                mooring_stiffness = _number(mooring["stiffness"],
                                            "mooring.stiffness", positive=True)
            elif mooring["kind"] == "moor_dyn":
                if set(mooring) != {"kind", "session", "point"}:
                    raise ValueError("MoorDyn needs a native session and spar-local point")
                moordyn, moordyn_point = mooring["session"], mooring["point"]
            else:
                raise ValueError("unsupported floating-joint mooring")
        imported_force = None
        auxiliary_files = ()
        if wave["type"] == "none":
            if set(wave) != {"type"} or "ramp_time" in sim:
                raise ValueError("no-wave floating joint has no wave or ramp settings")
            height, period, direction = 0.0, 8.0, 0.0
        elif wave["type"] == "elevationImport":
            steps = round(end_time / dt)
            if not np.isclose(steps * dt, end_time, rtol=0, atol=1e-10):
                raise ValueError("end_time must be an integer multiple of dt")
            time = np.arange(steps + 1) * dt
            elevation, imported_force, wave_path = _imported_elevation(
                wave, base, hydro[0], time, dt, ramp_time, rho, g,
            )
            auxiliary_files = (wave_path,)
            # The solver still requires nominal regular-wave arguments, but
            # the supplied excitation history drives this imported case.
            height, period, direction = 0.0, 8.0, 0.0
        elif wave["type"] in ("pm", "jonswap", "spectrumImport"):
            if wave["type"] == "spectrumImport":
                if (set(wave) != {"type", "file"}
                        or not isinstance(wave["file"], str) or not wave["file"]):
                    raise ValueError("floating-joint spectrumImport needs one MAT spectrum file")
                spectrum_path = (base / wave["file"]).resolve(strict=True)
                components = imported_spectrum_components(hydro[0], spectrum_path)
                auxiliary_files = (spectrum_path,)
            else:
                components, auxiliary_files = _irregular_components_from_case(
                    wave, hydro[0], base,
                )
            if (len(components.directions) != 1
                    or not np.isclose(components.directions[0], 0, atol=1e-12)):
                raise ValueError("floating-joint dynamics currently support one 0-degree wave heading")
            responses = tuple(synthesize_irregular_response(
                hydro[0], components, dt=dt, end_time=end_time,
                ramp_time=ramp_time, body_number=number, rho=rho, g=g,
                excitation_interpolation=wave.get("excitation_interpolation", "linear"),
            ) for number in (1, 2))
            elevation = responses[0].elevation
            imported_force = np.stack(
                [response.excitation_force for response in responses], axis=1,
            )
            height, period, direction = 0.0, 8.0, 0.0
        else:
            height = _number(wave.get("height"), "wave.height", nonnegative=True)
            period = _number(wave.get("period"), "wave.period", positive=True)
            if set(wave) - {"type", "height", "period", "direction"}:
                raise ValueError("regular waves use height, period, and direction")
            direction = _number(wave.get("direction", 0), "wave.direction")
            if direction != 0:
                raise ValueError("floating-joint dynamics currently support 0-degree waves")
        if wave["type"] == "regular" and "radiation_memory" in sim:
            raise ValueError("regular-wave floating-joint dynamics use constant radiation")
        # This joint path defaults to constant damping for regular waves;
        # broad-band and prescribed-elevation cases use radiation memory.
        radiation_method = sim.get(
            "radiation_method",
            "constant" if wave["type"] == "regular" else "convolution",
        )
        if wave["type"] == "regular" and radiation_method != "constant":
            raise ValueError("regular waves need constant radiation")
        if (wave["type"] in ("regularCIC", "pm", "jonswap", "spectrumImport",
                             "none", "elevationImport")
                and radiation_method not in ("convolution", "fir", "state_space")):
            raise ValueError("radiation memory needs convolution, FIR, or fitted radiation")
        if wave["type"] == "elevationImport" and radiation_method != "convolution":
            raise ValueError("imported elevation currently requires convolution radiation")
        if wave["type"] in ("regularCIC", "pm", "jonswap", "spectrumImport",
                            "none", "elevationImport"):
            radiation_memory = _number(
                sim.get("radiation_memory", 60), "simulation.radiation_memory",
                positive=True,
            )
        else:
            radiation_memory = None
        inertias = tuple(
            _number(body.get("pitch_inertia"), "body.pitch_inertia", positive=True)
            for body in bodies
        )
        solved = solve_rm3_regular(
            hydro[0], wave_height=height, wave_period=period,
            pitch_inertias=inertias, pto_damping=damping,
            pto_stiffness=stiffness, pto_equilibrium=equilibrium,
            pto_hard_stops=hard_stops,
            mooring_surge_stiffness=mooring_stiffness,
            moordyn=moordyn, moordyn_point=moordyn_point,
            b2b=b2b, radiation_memory=radiation_memory,
            radiation_method=radiation_method,
            added_mass_scheme=sim.get("added_mass_scheme", "implicit"),
            solver=sim.get("solver", "implicit"),
            no_wave=wave["type"] == "none",
            excitation_force=imported_force,
            initial_coordinate=initial_q, initial_speed=initial_v,
            joint_z=location[2],
            dt=dt, end_time=end_time, ramp_time=ramp_time, rho=rho, g=g,
        )
        if wave["type"] == "none":
            elevation = None
        elif wave["type"] not in ("elevationImport", "pm", "jonswap",
                                 "spectrumImport"):
            ramp = np.ones(len(solved.time))
            if ramp_time > 0:
                early = solved.time < ramp_time
                ramp[early] = (1 - np.cos(np.pi * solved.time[early] / ramp_time)) / 2
            elevation = height / 2 * ramp * np.cos(2 * np.pi * solved.time / period)
        extra_outputs = []
        if hard_stops is not None:
            extra_outputs.append(("pto_stop_force", solved.pto_stop_force))
        if mooring_stiffness:
            extra_outputs.extend((
                ("mooring_surge_position", solved.mooring_surge_position),
                ("mooring_surge_force", solved.mooring_surge_force),
            ))
        if moordyn is not None:
            extra_outputs.extend((
                ("moordyn_connection_position", solved.moordyn_connection_position),
                ("moordyn_connection_velocity", solved.moordyn_connection_velocity),
                ("moordyn_connection_force", solved.moordyn_connection_force),
            ))
        return CaseResponse(
            solved.time, solved.body_position, solved.body_velocity, hydro,
            pto_force=solved.pto_force, pto_label="pto_relative_heave_force",
            wave_elevation=elevation,
            auxiliary_files=auxiliary_files,
            extra_outputs=tuple(extra_outputs),
            coordinate_position=solved.coordinate_position,
            coordinate_velocity=solved.coordinate_velocity,
            pto_stroke=solved.pto_stroke,
            pto_velocity=solved.pto_velocity,
            pto_absorbed_power=solved.pto_dissipated_power,
        )

    raise ValueError(f"unsupported constraint layout: {kind}")

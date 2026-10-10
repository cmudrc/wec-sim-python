"""Case-driven entry point for the currently supported device layouts.

The case describes wave, body, constraint, and PTO properties rather than
naming RM3, OSWEC, or Sphere. Dynamic layouts feed the common
``GeneralizedDynamics`` engine through a validated adapter; a fixed Morison
layout evaluates stationary-body forces. Unsupported physics fails explicitly.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Mapping

import numpy as np
from scipy.io import loadmat

from .bodyClass import BodyClass
from .directDrive import integrate_direct_drive_heave
from .directLinearGenerator import DirectLinearGenerator
from .floatingOwc import FloatingOwcChamber, FloatingOwcTurbine
from .floatingOwcDynamics import solve_floating_owc
from .generalDynamics import BodyMotion, DynamicBody, GeneralizedDynamics
from .gbmFloating import (
    solve_floating_gbm_pm_orifice, solve_floating_gbm_regular,
)
from .hardStops import LinearHardStops
from .hingePitch import (
    solve_hinged_pitch_from_excitation, solve_hinged_pitch_regular,
)
from .irregularWave import (
    imported_full_directional_components, imported_spectrum_components,
    jonswap_equal_energy_components,
    pm_equal_energy_components,
    synthesize_full_directional_response, synthesize_irregular_response,
    synthesize_multiple_irregular_response,
)
from .linearCoordinates import build_coordinate_maps, initial_coordinate
from .linearHeave import solve_heave_free_decay
from .morison import (
    MorisonElement, no_wave_heave_morison_terms,
    regular_wave_axial_morison_terms, regular_wave_heave_morison_terms,
    solve_fixed_morison_irregular,
)
from .moorDyn import MoorDyn
from .nonlinearHydro import HeaveMeshHydro
from .orifice import OrificePTO
from .passiveYaw import (
    HeldPassiveYawExcitation, NearestHeadingExcitation,
    NearestSampledHeadingExcitation,
    PassiveYawExcitation, SampledPassiveYawExcitation,
)
from .ptoConnections import build_linear_ptos
from .rm3Regular import solve_rm3_regular
from .variableHydro import HeaveHydroState, integrate_variable_heave


AXES = ("surge", "sway", "heave", "roll", "pitch", "yaw")


@dataclass(frozen=True)
class CaseResponse:
    time: np.ndarray
    body_position: np.ndarray  # (time, body, six WEC-Sim coordinates)
    body_velocity: np.ndarray
    hydro_files: tuple[Path | None, ...]
    pto_force: np.ndarray | None = None
    pto_label: str | None = None
    wave_elevation: np.ndarray | None = None
    total_heave_force: np.ndarray | None = None
    auxiliary_files: tuple[Path, ...] = ()
    pto_generalized_force: np.ndarray | None = None
    extra_outputs: tuple[tuple[str, np.ndarray], ...] = ()
    coordinate_position: np.ndarray | None = None
    coordinate_velocity: np.ndarray | None = None
    pto_stroke: np.ndarray | None = None
    pto_velocity: np.ndarray | None = None
    pto_absorbed_power: np.ndarray | None = None


def _section(value, name, required, allowed):
    if not isinstance(value, Mapping):
        raise ValueError(f"{name} must be an object")
    missing = required - value.keys()
    extra = value.keys() - allowed
    if missing or extra:
        raise ValueError(f"{name}: missing {sorted(missing)}, unsupported {sorted(extra)}")
    return value


def _number(value, name, *, positive=False, nonnegative=False):
    if isinstance(value, bool):
        raise ValueError(f"{name} must be numeric")
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be numeric") from exc
    if (not np.isfinite(number) or (positive and number <= 0)
            or (nonnegative and number < 0)):
        raise ValueError(f"{name} is outside its supported range")
    return number


def _morison_elements(specs):
    if not isinstance(specs, list):
        raise ValueError("morison_elements must be a list")
    elements = []
    for spec in specs:
        spec = _section(
            spec, "morison element",
            {"point", "drag_coefficient", "added_mass_coefficient",
             "area", "volume"},
            {"point", "drag_coefficient", "added_mass_coefficient",
             "area", "volume", "phase_mode"},
        )
        elements.append(MorisonElement(
            tuple(spec["point"]), tuple(spec["drag_coefficient"]),
            tuple(spec["added_mass_coefficient"]), tuple(spec["area"]),
            spec["volume"], spec.get("phase_mode", "directional"),
        ))
    return tuple(elements)


def _hydro_file(body, base_dir):
    if isinstance(body, Mapping) and body.get("nonhydro", False):
        _section(body, "fixed nonhydrodynamic body",
                 {"nonhydro", "fixed", "center_gravity"},
                 {"nonhydro", "fixed", "center_gravity", "name", "mass", "inertia",
                  "volume", "morison_elements"})
        if body["nonhydro"] is not True or body["fixed"] is not True:
            raise ValueError("nonhydrodynamic body must be fixed in this layout")
        center = np.asarray(body["center_gravity"], dtype=float)
        if center.shape != (3,) or not np.isfinite(center).all():
            raise ValueError("fixed body center_gravity must be three finite coordinates")
        if "mass" in body:
            if body["mass"] != "equilibrium":
                _number(body["mass"], "fixed body mass", positive=True)
        if "inertia" in body:
            inertia = np.asarray(body["inertia"], dtype=float)
            if inertia.shape != (3,) or not np.isfinite(inertia).all() or np.any(inertia <= 0):
                raise ValueError("fixed body inertia must have three positive values")
        return None
    _section(body, "body", {"hydro_file"},
             {"hydro_file", "hydro_body", "mass", "pitch_inertia",
              "inertia", "coordinate_map", "name", "mean_drift", "fixed",
              "passive_yaw", "passive_yaw_threshold", "yaw_heading_bank",
              "geometry_file", "nonlinear_hydro",
              "drag_coefficient", "drag_area", "variable_hydro", "morison_elements"})
    raw = body["hydro_file"]
    if not isinstance(raw, str) or not raw:
        raise ValueError("body.hydro_file must be a file path")
    path = (base_dir / raw).expanduser().resolve(strict=True)
    if not path.is_file():
        raise ValueError(f"body.hydro_file is not a file: {path}")
    return path


def _body_number(body, expected):
    number = body.get("hydro_body", expected)
    if isinstance(number, bool) or number != expected:
        raise ValueError(f"body.hydro_body must be {expected} for this layout")


def _location(constraint, name="constraint.location"):
    location = np.asarray(constraint.get("location", [0, 0, 0]), dtype=float)
    if location.shape != (3,) or not np.isfinite(location).all():
        raise ValueError(f"{name} must have three finite coordinates")
    if not np.isclose(location[:2], 0, atol=1e-10).all():
        raise ValueError("supported joints must lie on the body x=y=0 axis")
    return location


def _pto(pto, kind, *, allow_location=False):
    _section(pto, "pto", {"kind", "damping"},
             {"kind", "damping", "stiffness", "equilibrium_position",
              "pretension"} | ({"location"} if allow_location else set()))
    if pto["kind"] != kind:
        raise ValueError(f"this layout requires a {kind} PTO")
    if "equilibrium_position" in pto and "pretension" in pto:
        raise ValueError("specify either PTO equilibrium_position or pretension")
    damping = _number(pto["damping"], "pto.damping", nonnegative=True)
    stiffness = _number(pto.get("stiffness", 0), "pto.stiffness", nonnegative=True)
    if "pretension" in pto:
        pretension = _number(pto["pretension"], "pto.pretension")
        if pretension and not stiffness:
            raise ValueError("nonzero PTO pretension needs positive stiffness")
        equilibrium = -pretension / stiffness if stiffness else 0.0
    else:
        equilibrium = _number(pto.get("equilibrium_position", 0),
                              "pto.equilibrium_position")
    if not np.isfinite(equilibrium):
        raise ValueError("PTO equilibrium position is outside its supported range")
    if equilibrium and not stiffness:
        raise ValueError("nonzero PTO equilibrium_position needs positive stiffness")
    return damping, stiffness, equilibrium


def _hard_stops(spec):
    required = {"lower_bound", "upper_bound", "lower_stiffness", "upper_stiffness"}
    optional = {"lower_damping", "upper_damping", "lower_transition_width",
                "upper_transition_width"}
    _section(spec, "pto.hard_stops", required, required | optional)
    values = {
        name: _number(value, f"pto.hard_stops.{name}",
                      positive=name.endswith("stiffness") or name.endswith("width"),
                      nonnegative=name.endswith("damping"))
        for name, value in spec.items()
    }
    return LinearHardStops(**values)


def _sampled_elevation(wave, base, time, ramp_time):
    """Load and ramp a sampled MATLAB elevation record."""
    _section(wave, "imported wave", {"type", "file"},
             {"type", "file", "variable", "direction", "reapply_force_ramp"})
    raw = wave["file"]
    if not isinstance(raw, (str, Path)) or not str(raw):
        raise ValueError("wave.file must be a MAT file path")
    path = (base / raw).resolve(strict=True)
    if not path.is_file():
        raise ValueError(f"wave.file is not a file: {path}")
    variable = wave.get("variable", "etaData")
    if not isinstance(variable, str) or not variable:
        raise ValueError("wave.variable must be a nonempty MAT variable name")
    direction = _number(wave.get("direction", 0), "wave.direction")
    if direction != 0:
        raise ValueError("imported elevation currently supports 0-degree waves")
    second_ramp = wave.get("reapply_force_ramp", False)
    if not isinstance(second_ramp, bool):
        raise ValueError("wave.reapply_force_ramp must be a boolean")
    mat = loadmat(path, variable_names=[variable])
    if variable not in mat:
        raise ValueError(f"wave.variable {variable!r} is absent from {path}")
    samples = np.asarray(mat[variable])
    if not np.issubdtype(samples.dtype, np.number) or not np.isrealobj(samples):
        raise ValueError("imported elevation must contain real numeric samples")
    samples = samples.astype(float, copy=False)
    if (samples.ndim != 2 or samples.shape[1] != 2 or samples.shape[0] < 2
            or not np.isfinite(samples).all()
            or not np.all(np.diff(samples[:, 0]) > 0)):
        raise ValueError("imported elevation must be a finite, increasing N-by-2 time/elevation array")
    if samples[0, 0] > time[0] + 1e-10 or samples[-1, 0] < time[-1] - 1e-10:
        raise ValueError("imported elevation must cover the full simulation time")
    ramp = np.ones_like(time)
    if ramp_time > 0:
        early = time < ramp_time
        ramp[early] = (1 - np.cos(np.pi * time[early] / ramp_time)) / 2
    elevation = np.interp(time, samples[:, 0], samples[:, 1]) * ramp
    return elevation, ramp, path


def _imported_elevation(wave, base, hydro_file, time, dt, ramp_time, rho, g):
    """Build RM3 body forces from one sampled MATLAB elevation record."""
    elevation, ramp, path = _sampled_elevation(wave, base, time, ramp_time)
    direction = _number(wave.get("direction", 0), "wave.direction")
    force = np.zeros((len(time), 2, 6))
    for number in (1, 2):
        body = BodyClass(str(hydro_file))
        body.bodyNumber = number
        body.bodyTotal = 2
        body.readH5file()
        body.hydroForce["userDefinedFe"] = np.zeros((len(time), 6))
        body.userDefinedExcitation(np.vstack((time, elevation)), dt, [direction], rho, g)
        force[:, number - 1] = body.hydroForce["userDefinedFe"]
    if wave.get("reapply_force_ramp", False):
        # The pinned MATLAB body block ramps force after waveClass ramped elevation.
        force *= ramp[:, None, None]
    return elevation, force, path


def _irregular_components_from_case(wave, hydro_file, base_dir):
    if "water_depth" in wave:
        raise ValueError("wave.water_depth override currently needs a fixed Morison body")
    if set(wave) - {"type", "height", "period", "directions", "spreading",
                     "seed", "phase_file", "phase_generator", "frequency_count",
                     "excitation_interpolation", "gamma", "frequency_range",
                     "discretization"}:
        raise ValueError("irregular waves use height, period, directions, and phase settings")
    if wave["type"] == "pm" and "gamma" in wave:
        raise ValueError("gamma applies only to JONSWAP waves")
    height = _number(wave.get("height"), "wave.height", positive=True)
    period = _number(wave.get("period"), "wave.period", positive=True)
    if "seed" in wave and "phase_file" in wave:
        raise ValueError("supply either wave.seed or wave.phase_file")
    auxiliary = ()
    if "phase_file" in wave:
        if not isinstance(wave["phase_file"], str) or not wave["phase_file"]:
            raise ValueError("wave.phase_file must be a file path")
        phase_path = (base_dir / wave["phase_file"]).resolve(strict=True)
        phase = np.loadtxt(phase_path, delimiter=",", ndmin=2)
        auxiliary = (phase_path,)
        seed = None
    else:
        phase = None
        seed = wave.get("seed", 7)
        if not isinstance(seed, int) or isinstance(seed, bool):
            raise ValueError("wave.seed must be an integer")
    builder = (jonswap_equal_energy_components if wave["type"] == "jonswap"
               else pm_equal_energy_components)
    options = {"gamma": wave["gamma"]} if "gamma" in wave else {}
    components = builder(
        hydro_file, significant_height=height, peak_period=period,
        directions=wave.get("directions", [0.0]),
        spreading=wave.get("spreading", [1.0]),
        count=wave.get("frequency_count"), seed=seed, phase=phase,
        phase_generator=wave.get("phase_generator", "numpy"),
        frequency_range=wave.get("frequency_range"),
        discretization=wave.get("discretization", "equal_energy"),
        **options,
    )
    return components, auxiliary


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
                    "radiation_memory", "radiation_method", "added_mass_scheme"})
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
            extra_outputs=tuple(extras),
        )

    if kind == "linear_subspace":
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
        if (len(bodies) not in (1, 2)
                or wave["type"] not in ("pm", "pm_multi", "regular",
                                        "spectrumImportFullDir")
                or b2b):
            raise ValueError("fixed-hinge pitch needs one flap, optional fixed base, and PM or regular waves")
        if len(bodies) == 2 and hydro[1] is not None:
            _body_number(bodies[1], 2)
            fixed_body = BodyClass(str(hydro[1]))
            fixed_body.bodyNumber = 2
            fixed_body.readH5file()
            fixed_center = np.asarray(fixed_body.cg, dtype=float).ravel()
        elif len(bodies) == 2:
            fixed_center = np.asarray(bodies[1]["center_gravity"], dtype=float)
        _body_number(bodies[0], 1)
        if set(bodies[0]) - {"hydro_file", "hydro_body", "mass", "pitch_inertia", "name"}:
            raise ValueError("fixed-hinge pitch uses mass and pitch_inertia")
        mass = _number(bodies[0].get("mass"), "body.mass", positive=True)
        inertia = _number(bodies[0].get("pitch_inertia"),
                          "body.pitch_inertia", positive=True)
        if set(constraint) - {"kind", "location"}:
            raise ValueError("fixed-hinge initial conditions are not yet supported")
        location = _location(constraint)
        pto_data = case.get("pto")
        damping, stiffness, equilibrium = _pto(
            pto_data, "pitch", allow_location=True,
        )
        if len(bodies) == 2 and "location" not in pto_data:
            raise ValueError("fixed nonhydrodynamic base needs pto.location")
        hinge = (_location({"location": pto_data["location"]}, "pto.location")
                 if "location" in pto_data else location)
        if wave["type"] == "regular":
            if sim.get("added_mass_scheme", "implicit") != "implicit":
                raise ValueError("regular fixed-hinge dynamics need implicit added mass")
            height = _number(wave.get("height"), "wave.height", positive=True)
            period = _number(wave.get("period"), "wave.period", positive=True)
            if set(wave) - {"type", "height", "period", "direction"}:
                raise ValueError("regular waves use height, period, and direction")
            direction = _number(wave.get("direction", 0), "wave.direction")
            if direction != 0:
                raise ValueError("fixed-hinge regular dynamics currently support 0-degree waves")
            if "radiation_memory" in sim:
                raise ValueError("regular fixed-hinge dynamics use constant-frequency radiation")
            solved = solve_hinged_pitch_regular(
                hydro[0], wave_height=height, wave_period=period,
                hinge_z=hinge[2], body_mass=mass, pitch_inertia=inertia,
                pto_damping=damping, pto_stiffness=stiffness,
                pto_equilibrium=equilibrium, dt=dt, end_time=end_time,
                ramp_time=ramp_time, rho=rho, g=g,
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
            if set(wave) - {"type", "height", "period", "directions",
                             "spreading", "seed", "phase_file", "phase_generator",
                             "frequency_count", "excitation_interpolation"}:
                raise ValueError("PM waves use height, period, directions, and phase settings")
            height = _number(wave.get("height"), "wave.height", positive=True)
            period = _number(wave.get("period"), "wave.period", positive=True)
            if "direction" in wave:
                raise ValueError("PM waves use directions and spreading arrays")
            if "seed" in wave and "phase_file" in wave:
                raise ValueError("supply either wave.seed or wave.phase_file")
            if "phase_file" in wave:
                if not isinstance(wave["phase_file"], str) or not wave["phase_file"]:
                    raise ValueError("wave.phase_file must be a file path")
                phase_path = (base / wave["phase_file"]).resolve(strict=True)
                phase = np.loadtxt(phase_path, delimiter=",", ndmin=2)
                seed = None
                auxiliary_files = [phase_path]
            else:
                phase = None
                seed = wave.get("seed", 7)
                if not isinstance(seed, int) or isinstance(seed, bool):
                    raise ValueError("wave.seed must be an integer")
                auxiliary_files = []
            count = wave.get("frequency_count", 500)
            components = pm_equal_energy_components(
                hydro[0], significant_height=height, peak_period=period,
                directions=wave.get("directions", [0, 30, 90]),
                spreading=wave.get("spreading", [0.1, 0.2, 0.7]),
                count=count, seed=seed, phase=phase,
                phase_generator=wave.get("phase_generator", "numpy"),
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
                "regular", "regularCIC", "pm", "jonswap", "none", "elevationImport"):
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
            height, period, direction = 0.0, 8.0, 0.0
        elif wave["type"] in ("pm", "jonswap"):
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
        radiation_method = sim.get(
            "radiation_method",
            "constant" if wave["type"] == "regular" else "convolution",
        )
        if wave["type"] == "regular" and radiation_method != "constant":
            raise ValueError("regular waves need constant radiation")
        if (wave["type"] in ("regularCIC", "pm", "jonswap", "none", "elevationImport")
                and radiation_method not in ("convolution", "fir")):
            raise ValueError("radiation memory needs convolution or FIR radiation")
        if wave["type"] == "elevationImport" and radiation_method != "convolution":
            raise ValueError("imported elevation currently requires convolution radiation")
        if wave["type"] in ("regularCIC", "pm", "jonswap", "none", "elevationImport"):
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
            no_wave=wave["type"] == "none",
            excitation_force=imported_force,
            initial_coordinate=initial_q, initial_speed=initial_v,
            joint_z=location[2],
            dt=dt, end_time=end_time, ramp_time=ramp_time, rho=rho, g=g,
        )
        if wave["type"] == "none":
            elevation = None
        elif wave["type"] not in ("elevationImport", "pm", "jonswap"):
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


def _current_settings(wave):
    if "current" not in wave:
        return 0.0, 0.0, "uniform", None
    current = _section(wave["current"], "wave.current",
                       {"speed", "direction", "profile"},
                       {"speed", "direction", "profile", "depth"})
    speed = _number(current["speed"], "wave.current.speed", nonnegative=True)
    direction = _number(current["direction"], "wave.current.direction")
    profile = current["profile"]
    if profile not in ("uniform", "power", "linear"):
        raise ValueError("wave.current.profile must be uniform, power, or linear")
    depth = (None if "depth" not in current else
             _number(current["depth"], "wave.current.depth", positive=True))
    if profile != "uniform" and depth is None:
        raise ValueError("depth-varying current needs wave.current.depth")
    return speed, direction, profile, depth


def _run_fixed_morison(case, sim, wave, constraint, bodies, hydro,
                       b2b, dt, end_time, ramp_time, rho, g, base_dir):
    """Evaluate the published type of stationary, no-HDF5 Morison device."""
    if (b2b or any(key in case for key in ("pto", "ptos", "mooring"))
            or set(constraint) != {"kind"}
            or set(sim) - {"dt", "end_time", "ramp_time", "rho", "g"}):
        raise ValueError("fixed Morison bodies use no joints, PTOs, or radiation")
    if (wave["type"] != "pm"
            or set(wave) - {"type", "height", "period", "directions",
                            "spreading", "seed", "phase_file", "phase_generator",
                            "frequency_count",
                            "frequency_range", "water_depth", "current"}
            or "frequency_range" not in wave or "water_depth" not in wave):
        raise ValueError("fixed Morison bodies need a PM sea, frequency range, and water depth")
    if "phase_file" in wave and "seed" in wave:
        raise ValueError("supply either wave.seed or wave.phase_file")
    auxiliary = ()
    if "phase_file" in wave:
        raw = wave["phase_file"]
        if not isinstance(raw, str) or not raw:
            raise ValueError("wave.phase_file must be a file path")
        phase_path = (base_dir / raw).resolve(strict=True)
        phase = np.loadtxt(phase_path, delimiter=",", ndmin=2)
        seed = None
        auxiliary = (phase_path,)
    else:
        phase = None
        seed = wave.get("seed", 7)
        if not isinstance(seed, int) or isinstance(seed, bool):
            raise ValueError("wave.seed must be an integer")
    depth = _number(wave["water_depth"], "wave.water_depth", positive=True)
    speed, current_direction, profile, current_depth = _current_settings(wave)
    components = pm_equal_energy_components(
        None,
        significant_height=_number(wave.get("height"), "wave.height", positive=True),
        peak_period=_number(wave.get("period"), "wave.period", positive=True),
        directions=wave.get("directions", [0.0]),
        spreading=wave.get("spreading", [1.0]),
        count=wave.get("frequency_count", 500),
        seed=seed, phase=phase, frequency_range=wave["frequency_range"],
        phase_generator=wave.get("phase_generator", "numpy"),
    )
    time = np.arange(round(end_time / dt) + 1) * dt
    position = np.zeros((len(time), len(bodies), 6))
    velocity = np.zeros_like(position)
    force_outputs = []
    elevation = None
    if any(path is not None for path in hydro):
        raise ValueError("fixed Morison bodies cannot use HDF5 hydrodynamics")
    for index, body in enumerate(bodies):
        if body.get("nonhydro") is not True or body.get("fixed") is not True:
            raise ValueError("fixed Morison bodies must be nonhydrodynamic and fixed")
        center = np.asarray(body["center_gravity"], dtype=float)
        position[:, index, :3] = center
        name = body.get("name", f"body{index + 1}")
        if not isinstance(name, str) or not name or any(
                entry[0] == f"morison_force_{name}" for entry in force_outputs):
            raise ValueError("fixed Morison body names must be nonempty and unique")
        elements = _morison_elements(body.get("morison_elements", []))
        if elements:
            solved = solve_fixed_morison_irregular(
                components, elements, center_gravity=center,
                water_depth=depth, dt=dt, end_time=end_time,
                ramp_time=ramp_time, rho=rho, g=g,
                current_speed=speed, current_direction=current_direction,
                current_profile=profile, current_depth=current_depth,
            )
            elevation = solved.wave_elevation
            force = solved.force
        else:
            force = np.zeros((len(time), 6))
        force_outputs.append((f"morison_force_{name}", force))
    if elevation is None:
        raise ValueError("fixed Morison device needs at least one element")
    return CaseResponse(
        time, position, velocity, hydro, wave_elevation=elevation,
        auxiliary_files=auxiliary, extra_outputs=tuple(force_outputs),
    )


def _run_linear_subspace(case, sim, wave, constraint, bodies, hydro,
                         b2b, dt, end_time, ramp_time, rho, g, base_dir):
    """Run mapped coordinates with supported regular, sampled, or no waves."""
    if set(constraint) - {"kind", "initial_coordinate", "initial_speed",
                           "coordinates"}:
        raise ValueError("linear_subspace uses coordinate maps, not joint locations")
    if wave["type"] not in ("regular", "regularCIC", "pm", "jonswap",
                            "spectrumImport", "elevationImport", "none"):
        raise ValueError("linear_subspace supports regular, regularCIC, PM, JONSWAP, spectrumImport, elevationImport, or no waves")
    if b2b and len(set(hydro)) != 1:
        raise ValueError("body-to-body hydrodynamics need one shared HDF5 file")
    if wave["type"] == "none" and set(wave) != {"type"}:
        raise ValueError("no-wave cases have no wave height or period")
    if (wave["type"] in ("pm", "jonswap", "spectrumImport", "elevationImport")
            and any(body.get("mean_drift", "none") != "none" for body in bodies)):
        raise ValueError("irregular linear-subspace mean-drift forcing is not supported")
    components = None
    auxiliary_files = []
    if wave["type"] in ("regular", "regularCIC"):
        allowed = {"type", "height", "period", "direction"}
        if wave["type"] == "regular":
            allowed.add("current")
        if set(wave) - allowed:
            raise ValueError("regular waves use height, period, and direction")
        height = _number(wave.get("height"), "wave.height", nonnegative=True)
        period = _number(wave.get("period"), "wave.period", positive=True)
        direction = _number(wave.get("direction", 0), "wave.direction")
        frequency = 2 * np.pi / period
        if wave["type"] == "regular" and "radiation_memory" in sim:
            raise ValueError("regular-wave linear dynamics use constant radiation")
    elif wave["type"] in ("pm", "jonswap"):
        height = _number(wave.get("height"), "wave.height", positive=True)
        period = _number(wave.get("period"), "wave.period", positive=True)
        components, phase_files = _irregular_components_from_case(
            wave, hydro[0], base_dir,
        )
        auxiliary_files.extend(phase_files)
    elif wave["type"] == "spectrumImport":
        if set(wave) != {"type", "file"} or not isinstance(wave["file"], str) or not wave["file"]:
            raise ValueError("spectrumImport needs one MAT spectrum file")
        spectrum_file = (base_dir / wave["file"]).expanduser().resolve(strict=True)
        components = imported_spectrum_components(hydro[0], spectrum_file)
        auxiliary_files.append(spectrum_file)
    elif wave["type"] == "elevationImport":
        pass  # The sampled record is loaded after the simulation grid is built.
    else:
        if "ramp_time" in sim:
            raise ValueError("ramp_time is inapplicable to no-wave dynamics")
    if wave["type"] != "regular":
        memory_time = _number(
            sim.get("radiation_memory", 15), "simulation.radiation_memory",
            positive=True,
        )
        memory_steps = round(memory_time / dt)
        if not np.isclose(memory_steps * dt, memory_time, atol=1e-10):
            raise ValueError("radiation_memory must be a multiple of dt")
        convolution_time = np.arange(memory_steps + 1) * dt
    time = np.arange(round(end_time / dt) + 1) * dt
    if not np.isclose(time[-1], end_time, atol=1e-10):
        raise ValueError("end_time must be an integer multiple of dt")
    imported_elevation = None
    imported_ramp = None
    if wave["type"] == "elevationImport":
        imported_elevation, imported_ramp, wave_path = _sampled_elevation(
            wave, base_dir, time, ramp_time,
        )
        auxiliary_files.append(wave_path)
    body_names = []
    for index, body_spec in enumerate(bodies, start=1):
        _body_number(body_spec, index)
        if "pitch_inertia" in body_spec:
            raise ValueError("linear_subspace uses body.inertia, not pitch_inertia")
        name = body_spec.get("name", f"body{index}")
        if not isinstance(name, str) or not name or name in body_names:
            raise ValueError("body names must be nonempty and unique")
        body_names.append(name)
    loaded_bodies = []
    centers = []
    for index, path in enumerate(hydro, start=1):
        body = BodyClass(str(path))
        body.bodyNumber = index
        body.bodyTotal = len(bodies)
        drift_option = bodies[index - 1].get("mean_drift", "none")
        drift_flags = {"none": 0, "control_surface": 1,
                       "momentum_conservation": 2}
        if not isinstance(drift_option, str) or drift_option not in drift_flags:
            raise ValueError(
                "body.mean_drift must be none, control_surface, or momentum_conservation"
            )
        if drift_option != "none" and wave["type"] not in ("regular", "regularCIC"):
            raise ValueError("mean drift requires regular incident waves")
        body.meanDriftForce = drift_flags[drift_option]
        body.readH5file()
        if drift_option != "none":
            drift_data = body.hydroData["hydro_coeffs"]["mean_drift"]
            if (drift_data.ndim != 3 or drift_data.shape[0] != 6
                    or not np.isfinite(drift_data).all()):
                raise ValueError(
                    f"body{index} HDF5 lacks finite {drift_option} mean-drift coefficients"
                )
        if int(np.asarray(body.dof).item()) != 6:
            raise ValueError("linear_subspace needs six-DOF hydrodynamic bodies")
        center = np.asarray(body.cg, dtype=float).ravel()
        if center.shape != (3,):
            raise ValueError("body center must have three coordinates")
        loaded_bodies.append(body)
        centers.append(center)
    maps, coordinate_names = build_coordinate_maps(
        constraint, bodies, body_names, centers,
    )
    n = maps[0].shape[1]
    moving_morison = tuple(
        index for index, spec in enumerate(bodies)
        if spec.get("morison_elements")
    )
    moving_elements = ()
    moving_terms = None
    if "current" in wave and not moving_morison:
        raise ValueError("regular-wave current needs a moving Morison element")
    if moving_morison:
        heave_map = np.zeros((6, 1))
        heave_map[2, 0] = 1
        three_dof_map = np.zeros((6, 3))
        three_dof_map[0, 0] = 1
        three_dof_map[2, 1] = 1
        three_dof_map[4, 2] = 1
        heave_only = n == 1 and np.array_equal(maps[0], heave_map)
        surge_heave_pitch = (n == 3 and np.array_equal(maps[0], three_dof_map)
                             and wave["type"] == "regular")
        if (len(bodies) != 1 or moving_morison != (0,)
                or not (heave_only or surge_heave_pitch)
                or wave["type"] not in ("none", "regular") or b2b
                or any(key in case for key in ("pto", "ptos"))
                or bodies[0].get("nonlinear_hydro") is not None
                or bodies[0].get("passive_yaw", False)
                or bodies[0].get("mean_drift", "none") != "none"):
            raise ValueError(
                "moving Morison needs one hydrodynamic body in pure heave "
                "or regular-wave surge/heave/pitch without a PTO"
            )
        if wave["type"] == "regular" and direction != 0:
            raise ValueError("moving Morison needs zero-heading regular waves")
        if (wave["type"] == "regular"
                and not np.allclose(centers[0][:2], 0, rtol=0, atol=1e-10)):
            raise ValueError("moving Morison needs a body centered at x=y=0")
        if "current" in wave and not surge_heave_pitch:
            raise ValueError("moving Morison current needs regular-wave surge/heave/pitch")
        current_settings = _current_settings(wave)
        moving_elements = _morison_elements(bodies[0]["morison_elements"])
        if wave["type"] == "regular":
            raw_depth = np.asarray(loaded_bodies[0].hydroData[
                "simulation_parameters"]["water_depth"]).item()
            depth = (np.inf if str(raw_depth).lower() == "infinite"
                     else float(raw_depth))

            if surge_heave_pitch:
                def moving_terms(at_time, coordinate, speed):
                    pose = np.r_[centers[0], np.zeros(3)] + maps[0] @ coordinate
                    velocity = maps[0] @ speed
                    return regular_wave_axial_morison_terms(
                        moving_elements, position=pose, velocity=velocity,
                        time=at_time, wave_height=height, wave_period=period,
                        ramp_time=ramp_time, water_depth=depth, rho=rho, g=g,
                        current_speed=current_settings[0],
                        current_direction=current_settings[1],
                        current_profile=current_settings[2],
                        current_depth=current_settings[3],
                    )
            else:
                def moving_terms(at_time, coordinate, speed):
                    applied, added = regular_wave_heave_morison_terms(
                        moving_elements, center_z=centers[0][2],
                        heave=coordinate[0], speed=speed[0], time=at_time,
                        wave_height=height, wave_period=period,
                        ramp_time=ramp_time, water_depth=depth, rho=rho, g=g,
                    )
                    force = np.zeros(6)
                    force[2] = applied
                    matrix = np.zeros((6, 6))
                    matrix[2, 2] = added
                    return force, matrix
        else:
            def moving_terms(at_time, coordinate, speed):
                applied, added = no_wave_heave_morison_terms(
                    moving_elements, center_z=centers[0][2],
                    heave=coordinate[0], speed=speed[0], rho=rho,
                )
                force = np.zeros(6)
                force[2] = applied
                matrix = np.zeros((6, 6))
                matrix[2, 2] = added
                return force, matrix
        moving_terms(0, np.zeros(n), np.zeros(n))
    passive_indices = [
        index for index, spec in enumerate(bodies)
        if spec.get("passive_yaw", False) is True
    ]
    if any(not isinstance(spec.get("passive_yaw", False), bool)
           for spec in bodies):
        raise ValueError("body.passive_yaw must be a boolean")
    yaw_thresholds = [
        _number(spec.get("passive_yaw_threshold", 0),
                "body.passive_yaw_threshold", nonnegative=True)
        for spec in bodies
    ]
    if any(threshold and (not bodies[index].get("passive_yaw", False)
                          or wave["type"] != "pm")
           for index, threshold in enumerate(yaw_thresholds)):
        raise ValueError("positive passive_yaw_threshold needs PM passive yaw")
    yaw_banks = [spec.get("yaw_heading_bank") for spec in bodies]
    if any(bank is not None and (not bodies[index].get("passive_yaw", False)
                                 or yaw_thresholds[index]
                                 or wave["type"] not in ("regular", "pm"))
           for index, bank in enumerate(yaw_banks)):
        raise ValueError("yaw_heading_bank needs regular or PM passive yaw without a threshold")
    if passive_indices:
        yaw_map = np.zeros((6, 1))
        yaw_map[5, 0] = 1
        if (len(passive_indices) != 1 or n != 1 or b2b
                or wave["type"] not in ("regular", "pm")
                or not np.array_equal(maps[passive_indices[0]], yaw_map)
                or any(np.any(mapping) for index, mapping in enumerate(maps)
                       if index != passive_indices[0])
                or bodies[passive_indices[0]].get("mean_drift", "none") != "none"):
            raise ValueError(
                "passive yaw currently needs one pure-yaw body, stationary others, "
                "regular or PM waves, and independent radiation"
            )
    nonlinear_indices = [
        index for index, spec in enumerate(bodies)
        if spec.get("nonlinear_hydro") is not None
    ]
    if nonlinear_indices:
        heave_map = np.zeros((6, 1))
        heave_map[2, 0] = 1
        if (len(bodies) != 1 or nonlinear_indices != [0] or n != 1
                or not np.array_equal(maps[0], heave_map)
                or wave["type"] not in ("regular", "regularCIC")
                or direction != 0 or b2b
                or bodies[0]["nonlinear_hydro"] != "instantaneous"
                or bodies[0].get("mean_drift", "none") != "none"
                or bodies[0].get("passive_yaw", False)):
            raise ValueError(
                "instantaneous nonlinear hydro currently needs one pure-heave "
                "body and a zero-direction regular or regularCIC wave"
            )
    for spec in bodies:
        if spec.get("nonlinear_hydro") is None and any(
            key in spec for key in ("geometry_file", "drag_coefficient", "drag_area")
        ):
            raise ValueError("mesh geometry and drag settings require nonlinear_hydro")
    initial_q = initial_coordinate(
        constraint.get("initial_coordinate", [0] * n), coordinate_names,
        "constraint.initial_coordinate",
    )
    initial_v = initial_coordinate(
        constraint.get("initial_speed", [0] * n), coordinate_names,
        "constraint.initial_speed",
    )
    if "pto" in case and "ptos" in case:
        raise ValueError("use either a PTO matrix or PTO connections")
    if "pto" in case:
        pto = _section(case["pto"], "pto", {"kind"},
                       {"kind", "stiffness_matrix", "damping_matrix",
                        "equilibrium_coordinate"})
        if pto["kind"] != "linear":
            raise ValueError("linear_subspace requires a linear PTO matrix")
        stiffness = np.asarray(pto.get("stiffness_matrix", np.zeros((n, n))), dtype=float)
        damping = np.asarray(pto.get("damping_matrix", np.zeros((n, n))), dtype=float)
        if (stiffness.shape != (n, n) or damping.shape != (n, n)
                or not np.isfinite(stiffness).all() or not np.isfinite(damping).all()):
            raise ValueError("PTO matrices must be finite N-by-N matrices")
        equilibrium = np.asarray(pto.get("equilibrium_coordinate", [0] * n),
                                 dtype=float)
        if equilibrium.shape != (n,) or not np.isfinite(equilibrium).all():
            raise ValueError("pto.equilibrium_coordinate must have N finite values")
        if np.any(equilibrium) and not np.any(stiffness @ equilibrium):
            raise ValueError("PTO equilibrium offset must produce a spring force")
    else:
        stiffness = np.zeros((n, n))
        damping = np.zeros((n, n))
        equilibrium = np.zeros(n)
    pto_bias = np.zeros(n)
    connections = ()

    dynamic_bodies = []
    nonlinear_models = []
    passive_model = None
    pm_elevation = None
    for index, (body_spec, body, mapping) in enumerate(
            zip(bodies, loaded_bodies, maps), start=1):
        mass_setting = body_spec.get("mass", "equilibrium")
        if mass_setting != "equilibrium":
            mass_setting = _number(mass_setting, "body.mass", positive=True)
        mesh_model = None
        if index - 1 in nonlinear_indices:
            geometry_file = body_spec.get("geometry_file")
            if not isinstance(geometry_file, str) or not geometry_file:
                raise ValueError("instantaneous nonlinear hydro needs geometry_file")
            geometry_path = (base_dir / geometry_file).expanduser().resolve(strict=True)
            cd = _number(body_spec.get("drag_coefficient", 0),
                         "body.drag_coefficient", nonnegative=True)
            drag_area = _number(body_spec.get("drag_area", 0),
                                "body.drag_area", nonnegative=True)
            depth = np.asarray(
                body.hydroData["simulation_parameters"]["water_depth"]
            ).ravel()
            if depth.size != 1 or not np.isfinite(depth[0]) or depth[0] <= 0:
                raise ValueError("nonlinear hydro needs positive HDF5 water depth")
            if not np.allclose(centers[index - 1][:2], 0, rtol=0, atol=1e-10):
                raise ValueError("heave mesh hydro currently needs CG at x=y=0")
            mesh_model = HeaveMeshHydro.from_stl(
                geometry_path, center_z=centers[index - 1][2], rho=rho,
                gravity=g, depth=float(depth[0]), period=period,
                height=height, ramp_time=ramp_time,
                mass=None if mass_setting == "equilibrium" else mass_setting,
                drag_coefficient=cd, drag_area=drag_area,
            )
            if mass_setting == "equilibrium":
                mass_setting = mesh_model.mass
            auxiliary_files.append(geometry_path)
        nonlinear_models.append(mesh_model)
        body.mass = mass_setting
        inertia = np.asarray(body_spec.get("inertia", [0, 0, 0]), dtype=float)
        if (inertia.shape != (3,) or not np.isfinite(inertia).all()
                or np.any(inertia < 0)):
            raise ValueError("body.inertia must have three nonnegative entries")
        body.hydroStiffness = np.zeros((6, 6))
        body.viscDrag = {
            "Drag": np.zeros((6, 6)), "cd": np.zeros(6),
            "characteristicArea": np.zeros(6),
        }
        body.linearDamping = np.zeros((6, 6))
        wave_amp = np.vstack((
            time, imported_elevation if imported_elevation is not None
            else np.zeros_like(time),
        ))
        if wave["type"] == "regular":
            body.hydroForcePre(
                frequency, [direction], 1, np.array([0.0]), [], dt, rho, g,
                "regular", wave_amp, index, len(bodies), 0, 0, int(b2b),
            )
        else:
            irf_time = body.hydroData["hydro_coeffs"]["radiation_damping"][
                "impulse_response_fun"]["t"]
            if memory_time > np.max(irf_time) + 1e-10:
                raise ValueError("radiation_memory exceeds the HDF5 kernel")
            regular_memory = wave["type"] == "regularCIC"
            irregular = wave["type"] in ("pm", "jonswap", "spectrumImport")
            body.hydroForcePre(
                (frequency if regular_memory else
                 components.omega if irregular else []),
                ([direction] if regular_memory else
                 components.directions if irregular else [0]),
                len(convolution_time), convolution_time,
                len(components.omega) if irregular else [],
                dt, rho, g, ("regularCIC" if regular_memory else
                             "spectrumImport" if wave["type"] == "spectrumImport" else
                             "elevationImport" if wave["type"] == "elevationImport" else
                             "irregular" if irregular else "noWaveCIC"),
                wave_amp,
                index, len(bodies), 0, 0, int(b2b),
            )
        physical_mass = float(np.asarray(body.mass).item())
        rigid_mass = np.diag([physical_mass] * 3 + inertia.tolist())
        hydro_force = body.hydroForce
        if b2b:
            added_mass = tuple(
                np.asarray(hydro_force["fAddedMass"])[:, 6*j:6*(j+1)]
                for j in range(len(bodies))
            )
            if wave["type"] == "regular":
                radiation_damping = tuple(
                    np.asarray(hydro_force["fDamping"])[:, 6*j:6*(j+1)]
                    for j in range(len(bodies))
                )
        else:
            added_mass = [np.zeros((6, 6)) for _ in bodies]
            added_mass[index - 1] = np.asarray(hydro_force["fAddedMass"])
            if wave["type"] == "regular":
                radiation_damping = [np.zeros((6, 6)) for _ in bodies]
                radiation_damping[index - 1] = np.asarray(hydro_force["fDamping"])
        if wave["type"] != "regular":
            radiation_damping = tuple(np.zeros((6, 6)) for _ in bodies)
            kernel = np.asarray(hydro_force["irkb"])
            if not b2b and len(bodies) > 1:
                independent_kernel = np.zeros((len(kernel), 6, 6 * len(bodies)))
                independent_kernel[:, :, 6 * (index - 1):6 * index] = kernel
                kernel = independent_kernel
        else:
            kernel = None
        center = centers[index - 1]

        def motion(q, v, *, mapping=mapping):
            return BodyMotion(mapping @ q, mapping, np.zeros(6))

        if wave["type"] in ("regular", "regularCIC"):
            real = np.asarray(hydro_force["fExt"]["re"])
            imaginary = np.asarray(hydro_force["fExt"]["im"])
            drift = np.asarray(hydro_force["fExt"]["md"])

            def excitation(at_time, *, real=real, imaginary=imaginary,
                           drift=drift):
                ramp = (1.0 if ramp_time == 0 or at_time >= ramp_time
                        else (1 - np.cos(np.pi * at_time / ramp_time)) / 2)
                return (height / 2 * ramp * (
                    real * np.cos(frequency * at_time)
                    - imaginary * np.sin(frequency * at_time)
                ) + (height / 2)**2 * ramp * drift)
        elif wave["type"] == "elevationImport":
            sampled_force = np.asarray(hydro_force["userDefinedFe"])
            if wave.get("reapply_force_ramp", False):
                sampled_force = sampled_force * imported_ramp[:, None]

            def excitation(at_time, *, sampled_force=sampled_force):
                if len(time) == 1:
                    return sampled_force[0]
                sample = min(max(at_time / dt, 0.0), float(len(time) - 1))
                left = min(int(sample), len(time) - 2)
                fraction = sample - left
                return ((1 - fraction) * sampled_force[left]
                        + fraction * sampled_force[left + 1])
        elif wave["type"] in ("pm", "jonswap", "spectrumImport"):
            if body_spec.get("passive_yaw", False):
                def excitation(at_time):
                    return np.zeros(6)
            else:
                incident = synthesize_irregular_response(
                    hydro[index - 1], components, dt=dt, end_time=end_time,
                    ramp_time=ramp_time, body_number=index, rho=rho, g=g,
                    excitation_interpolation=wave.get(
                        "excitation_interpolation", "linear",
                    ),
                )
                if pm_elevation is None:
                    pm_elevation = incident.elevation
                sampled_force = incident.excitation_force

                def excitation(at_time, *, sampled_force=sampled_force):
                    sample = round(at_time / dt)
                    return sampled_force[sample]
        else:
            def excitation(at_time):
                return np.zeros(6)

        state_excitation = None
        if body_spec.get("passive_yaw", False):
            if wave["type"] == "pm":
                passive_model = SampledPassiveYawExcitation.from_hydro_data(
                    body.hydroData, components, dt=dt, end_time=end_time,
                    ramp_time=ramp_time, rho=rho, g=g,
                )
                pm_elevation = passive_model.elevation
                if yaw_banks[index - 1] is not None:
                    passive_model = NearestSampledHeadingExcitation(
                        passive_model, yaw_banks[index - 1],
                    )
                if yaw_thresholds[index - 1]:
                    passive_model = HeldPassiveYawExcitation(
                        passive_model, yaw_thresholds[index - 1],
                    )
            else:
                passive_model = PassiveYawExcitation.from_hydro_data(
                    body.hydroData, omega=frequency,
                    incident_direction=direction, amplitude=height / 2,
                    ramp_time=ramp_time, rho=rho, g=g,
                    spline_frequency=yaw_banks[index - 1] is not None,
                )
                if yaw_banks[index - 1] is not None:
                    passive_model = NearestHeadingExcitation(
                        passive_model, yaw_banks[index - 1],
                    )

            if isinstance(passive_model, HeldPassiveYawExcitation):
                state_excitation = passive_model
            else:
                def state_excitation(at_time, coordinate, speed, *,
                                     model=passive_model):
                    return model.force(at_time, coordinate[0])
        if mesh_model is not None:
            def state_excitation(at_time, coordinate, speed, *,
                                 model=mesh_model, linear=excitation):
                buoyancy, fk, drag = model.forces(
                    at_time, coordinate[0], speed[0]
                )
                result = linear(at_time).copy()
                result[2] += buoyancy + fk + drag
                return result

        state_inertia = None
        if moving_elements:
            def state_inertia(at_time, coordinate, speed, *, terms=moving_terms):
                return terms(at_time, coordinate, speed)

        dynamic_bodies.append(DynamicBody(
            rigid_mass=rigid_mass,
            added_mass=tuple(added_mass),
            damping=tuple(radiation_damping),
            restoring=(np.zeros((6, 6)) if mesh_model is not None else
                       np.asarray(hydro_force["linearHydroRestCoef"])),
            static_force=(np.zeros(6) if mesh_model is not None else
                          np.array([
                              0, 0, (rho * float(np.asarray(body.dispVol).item())
                                     - physical_mass) * g, 0, 0, 0,
                          ])),
            reference_position=np.r_[center, np.zeros(3)],
            motion=motion,
            excitation=excitation,
            radiation_kernel=kernel,
            state_excitation=state_excitation,
            state_inertia=state_inertia,
        ))
    if "ptos" in case:
        connections, stiffness, damping, pto_bias = build_linear_ptos(
            case["ptos"], maps, centers, body_names, coordinate_names,
            allow_linear_generator=True,
        )
    generator_specs = {
        spec["name"]: spec["linear_generator"]
        for spec in case.get("ptos", []) if "linear_generator" in spec
    }
    if generator_specs and wave["type"] != "regular":
        raise ValueError("linear generator currently needs regular-wave constant radiation")
    linear_generators = []
    for connection in connections:
        if connection.name not in generator_specs:
            continue
        name = f"pto.{connection.name}.linear_generator"
        spec = _section(
            generator_specs[connection.name], name,
            {"stator_resistance", "friction", "pole_pitch", "magnet_flux",
             "inductance", "load_resistance"},
            {"stator_resistance", "friction", "pole_pitch", "magnet_flux",
             "inductance", "load_resistance", "initial_angle",
             "initial_flux_d", "initial_flux_q"},
        )
        settings = {
            key: _number(value, f"{name}.{key}")
            for key, value in spec.items() if value is not None
        }
        linear_generators.append((connection, DirectLinearGenerator(**settings)))
    linear_generators = tuple(linear_generators)
    controlled_connections = tuple(
        connection for connection in connections if connection.control is not None
    )
    system = GeneralizedDynamics(
        tuple(dynamic_bodies), n, pto_stiffness=stiffness,
        pto_damping=damping, pto_equilibrium=equilibrium,
        pto_bias=pto_bias,
        controlled_ptos=controlled_connections,
        linear_generators=linear_generators,
    )
    response = system.integrate(
        dt=dt, end_time=end_time,
        initial_coordinate=initial_q, initial_speed=initial_v,
    )
    if wave["type"] in ("regular", "regularCIC"):
        ramp = np.ones(len(response.time))
        if ramp_time > 0:
            early = response.time < ramp_time
            ramp[early] = (1 - np.cos(np.pi * response.time[early] / ramp_time)) / 2
        elevation = height / 2 * ramp * np.cos(frequency * response.time)
    else:
        elevation = (imported_elevation if imported_elevation is not None
                     else pm_elevation)
    coordinate_outputs = tuple(
        output for index, name in enumerate(coordinate_names)
        for output in (
            (f"coordinate_{name}_position", response.coordinate[:, index]),
            (f"coordinate_{name}_velocity", response.speed[:, index]),
        )
    ) if "coordinates" in constraint else ()
    controlled_forces = (
        {connection.name: response.controlled_pto_force[:, index]
         for index, connection in enumerate(controlled_connections)}
        if controlled_connections else {}
    )
    generator_forces = (
        {connection.name: response.linear_generator_force[:, index]
         for index, (connection, _) in enumerate(linear_generators)}
        if linear_generators else {}
    )
    pto_outputs = tuple(
        output for connection in connections
        for output in connection.outputs(
            response.coordinate, response.speed,
            force_override=(controlled_forces | generator_forces).get(connection.name),
        )
    )
    generator_outputs = []
    for index, (connection, generator) in enumerate(linear_generators):
        states = response.linear_generator_state[:, index]
        stroke_speed = response.speed @ connection.stroke_jacobian
        electrical_power = np.empty(len(response.time))
        phase_current = np.empty((len(response.time), 3))
        phase_voltage = np.empty_like(phase_current)
        for step, (speed, state) in enumerate(zip(stroke_speed, states)):
            signals = generator.signals(float(speed), state)
            electrical_power[step] = signals.electrical_power
            phase_current[step] = signals.phase_current
            phase_voltage[step] = signals.phase_voltage
        prefix = f"pto_{connection.name}_generator_"
        generator_outputs.extend((
            (prefix + "flux_d", states[:, 0]),
            (prefix + "flux_q", states[:, 1]),
            (prefix + "angle", states[:, 2]),
            (prefix + "friction_force", generator.friction * stroke_speed),
            (prefix + "electrical_power", electrical_power),
            (prefix + "phase_current", phase_current),
            (prefix + "phase_voltage", phase_voltage),
        ))
    generalized_pto = (-(response.coordinate - equilibrium) @ stiffness.T
                       - response.speed @ damping.T + pto_bias)
    for connection in controlled_connections:
        generalized_pto += np.outer(
            controlled_forces[connection.name], connection.stroke_jacobian,
        )
    for connection, _ in linear_generators:
        generalized_pto += np.outer(
            generator_forces[connection.name], connection.stroke_jacobian,
        )
    drift_outputs = tuple(
        output
        for index, (spec, body, dynamic) in enumerate(
            zip(bodies, loaded_bodies, dynamic_bodies), start=1,
        )
        if spec.get("mean_drift", "none") != "none"
        for output in (
            (f"body{index}_mean_drift_force",
             (height / 2)**2 * ramp[:, None]
             * np.asarray(body.hydroForce["fExt"]["md"])[None, :]),
            (f"body{index}_excitation_force",
             np.stack([dynamic.excitation(t) for t in response.time])),
        )
    ) if wave["type"] in ("regular", "regularCIC") else ()
    imported_outputs = tuple(
        (f"body{index}_excitation_force",
         (np.asarray(body.hydroForce["userDefinedFe"])
          * (imported_ramp[:, None] if wave.get("reapply_force_ramp", False)
             else 1)))
        for index, body in enumerate(loaded_bodies, start=1)
    ) if wave["type"] == "elevationImport" else ()
    passive_outputs = ()
    if passive_model is not None:
        body_index = passive_indices[0] + 1
        if isinstance(passive_model, HeldPassiveYawExcitation):
            excitation_history = np.asarray(passive_model.force_history)
        else:
            excitation_history = np.stack([
                passive_model.force(t, angle)
                for t, angle in zip(response.time, response.coordinate[:, 0])
            ])
        passive_outputs = ((
            f"body{body_index}_excitation_force",
            excitation_history,
        ),)
    nonlinear_outputs = ()
    if nonlinear_indices:
        model = nonlinear_models[0]
        components = np.array([
            model.forces(t, q[0], v[0])
            for t, q, v in zip(response.time, response.coordinate, response.speed)
        ])
        nonlinear_outputs = (
            ("body1_buoyancy_minus_weight", components[:, 0]),
            ("body1_nonlinear_fk_correction", components[:, 1]),
            ("body1_quadratic_drag_force", components[:, 2]),
        )
    morison_outputs = ()
    if moving_elements:
        force = np.stack([
            applied - added @ (maps[0] @ a)
            for (applied, added), a in zip(
                (moving_terms(t, q, v)
                 for t, q, v in zip(response.time, response.coordinate,
                                    response.speed)),
                response.acceleration,
            )
        ])
        morison_outputs = ((f"morison_force_{body_names[0]}", force),)
    return CaseResponse(
        response.time, response.body_position, response.body_velocity,
        hydro, wave_elevation=elevation,
        auxiliary_files=tuple(auxiliary_files),
        pto_generalized_force=(
            generalized_pto if "pto" in case or "ptos" in case else None
        ),
        extra_outputs=(coordinate_outputs + pto_outputs + tuple(generator_outputs) + drift_outputs
                       + imported_outputs
                       + passive_outputs + nonlinear_outputs + morison_outputs),
    )

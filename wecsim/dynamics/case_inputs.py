"""Validated case input and wave helpers."""

from pathlib import Path
from typing import Mapping

import numpy as np
from scipy.io import loadmat

from wecsim.dynamics.hardStops import LinearHardStops
from wecsim.environment.irregularWave import (
    jonswap_equal_energy_components,
    pm_equal_energy_components,
)
from wecsim.hydrodynamics.bodyClass import BodyClass
from wecsim.hydrodynamics.morison import MorisonElement


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
             "area", "volume", "phase_mode", "axis"},
        )
        elements.append(MorisonElement(
            tuple(spec["point"]), tuple(spec["drag_coefficient"]),
            tuple(spec["added_mass_coefficient"]), tuple(spec["area"]),
            spec["volume"], spec.get("phase_mode", "directional"),
            tuple(spec["axis"]) if "axis" in spec else None,
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
            # This body never accelerates. Zero inertia is valid and is the
            # public fixed_body default; negative inertia is still invalid.
            if inertia.shape != (3,) or not np.isfinite(inertia).all() or np.any(inertia < 0):
                raise ValueError("fixed body inertia must have three nonnegative values")
        return None
    _section(body, "body", {"hydro_file"},
             {"hydro_file", "hydro_body", "mass", "pitch_inertia",
              "inertia", "coordinate_map", "name", "mean_drift", "fixed",
              "passive_yaw", "passive_yaw_threshold", "yaw_heading_bank",
              "geometry_file", "nonlinear_hydro",
              "drag_coefficient", "drag_area", "variable_hydro", "morison_elements",
              "hydrostatic_stiffness", "linear_damping", "quadratic_drag",
              "quadratic_drag_coefficient",
              "quadratic_drag_area"})
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
    # The MAT record has its own sample times. Interpolate onto the solver
    # grid before applying the same startup ramp used by generated waves.
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
    # A phase file replays a specific MATLAB realization; a seed asks the
    # selected phase generator to create a reproducible realization instead.
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

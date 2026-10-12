"""Fixed-body Morison case adapter."""

import numpy as np

from wecsim.dynamics.case_inputs import _current_settings, _morison_elements, _number
from wecsim.dynamics.case_response import CaseResponse
from wecsim.environment.irregularWave import pm_equal_energy_components
from wecsim.hydrodynamics.morison import solve_fixed_morison_irregular


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

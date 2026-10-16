"""Reject unsupported physics settings before they become misleading motion."""

from pathlib import Path

import pytest

from wecsim.dynamics.caseDynamics import run_case


HYDRO = str(Path(__file__).resolve().parents[1] / "examples/data/sphere.h5")


def _free_decay():
    return {
        "simulation": {"dt": .1, "end_time": 0},
        "wave": {"type": "none"},
        "bodies": [{"hydro_file": HYDRO}],
        "constraint": {"kind": "heave"},
    }


@pytest.mark.parametrize("section,key,value,reason", [
    ("wave", "current", {"speed": 1}, "wave.current requires"),
    ("simulation", "radiation_method", "state_space",
     "currently applies to floating_joint"),
    ("simulation", "solver", "wecsim_ode4",
     "needs linear_subspace or a nonlinear regular fixed_hinge"),
    ("simulation", "added_mass_scheme", "implicit", "currently applies to"),
    ("body", "fixed", True, "fixed hydrodynamic body"),
    ("body", "mean_drift", "control_surface", "requires linear_subspace"),
    ("case", "ptos", [], "require linear_subspace"),
    ("case", "mooring", {}, "requires a floating_joint or floating_owc"),
    ("case", "body_to_body", "yes", "must be a boolean"),
])
def test_free_decay_rejects_settings_from_other_physics(section, key, value,
                                                        reason):
    case = _free_decay()
    target = case["bodies"][0] if section == "body" else (
        case if section == "case" else case[section]
    )
    target[key] = value
    with pytest.raises(ValueError, match=reason):
        run_case(case)


def test_case_needs_at_least_one_body():
    case = _free_decay()
    case["bodies"] = []
    with pytest.raises(ValueError, match="bodies must be a nonempty list"):
        run_case(case)


def _mapped_heave():
    return {
        "simulation": {"dt": .1, "end_time": 0},
        "wave": {"type": "regular", "height": 0, "period": 8},
        "bodies": [{"hydro_file": HYDRO,
                    "coordinate_map": [[0], [0], [1], [0], [0], [0]]}],
        "constraint": {"kind": "linear_subspace"},
    }


@pytest.mark.parametrize("setting,reason", [
    ("no_wave_height", "no-wave cases have no wave height"),
    ("regular_memory", "regular-wave linear dynamics use constant radiation"),
    ("irregular_mean_drift", "irregular linear-subspace mean-drift"),
    ("current_without_morison", "current needs a moving Morison element"),
])
def test_mapped_motion_rejects_unsupported_force_combinations(setting, reason):
    case = _mapped_heave()
    if setting == "no_wave_height":
        case["wave"] = {"type": "none", "height": 1}
    elif setting == "regular_memory":
        case["simulation"]["radiation_memory"] = 1
    elif setting == "irregular_mean_drift":
        case["wave"]["type"] = "pm"
        case["bodies"][0]["mean_drift"] = "control_surface"
    else:
        case["wave"]["current"] = {"speed": 1}
    with pytest.raises(ValueError, match=reason):
        run_case(case)

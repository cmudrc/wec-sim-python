"""Keep the JSON simulation contract explicit for future validator changes."""

from copy import deepcopy
from pathlib import Path

import pytest

from wecsim.dynamics.caseDynamics import run_case


HYDRO = str(Path(__file__).resolve().parents[1] / "examples/data/sphere.h5")


def _regular_case():
    return {
        "simulation": {"dt": .1, "end_time": 0},
        "wave": {"type": "regular", "height": 0, "period": 8},
        "bodies": [{"hydro_file": HYDRO,
                    "coordinate_map": [[0], [0], [1], [0], [0], [0]]}],
        "constraint": {"kind": "linear_subspace"},
    }


def test_simulation_accepts_numeric_strings_without_mutating_case():
    case = _regular_case()
    case["simulation"] = {
        "dt": "0.1", "end_time": "0", "ramp_time": "0",
        "rho": "1000", "g": "9.81",
    }
    original = deepcopy(case)

    response = run_case(case)

    assert response.time.tolist() == [0]
    assert case == original


@pytest.mark.parametrize("update, message", [
    ({"dt": True}, "simulation.dt must be numeric"),
    ({"end_time": False}, "simulation.end_time must be numeric"),
    ({"dt": "bad"}, "simulation.dt must be numeric"),
    ({"dt": "0"}, "simulation.dt is outside its supported range"),
    ({"end_time": "-1"}, "simulation.end_time is outside its supported range"),
    ({"ramp_time": "-1"}, "simulation.ramp_time is outside its supported range"),
    ({"rho": False}, "simulation.rho must be numeric"),
    ({"rho": "NaN"}, "simulation.rho is outside its supported range"),
    ({"g": "Infinity"}, "simulation.g is outside its supported range"),
    ({"unknown": 1}, "simulation: missing [], unsupported ['unknown']"),
    ({"dt": True, "end_time": "bad"}, "simulation.dt must be numeric"),
])
def test_simulation_rejects_invalid_values_and_unknown_keys(update, message):
    case = _regular_case()
    case["simulation"].update(update)
    with pytest.raises(ValueError) as exc:
        run_case(case)
    assert str(exc.value) == message


def test_simulation_requires_object_and_required_keys():
    case = _regular_case()
    case["simulation"] = []
    with pytest.raises(ValueError) as exc:
        run_case(case)
    assert str(exc.value) == "simulation must be an object"

    case["simulation"] = {"dt": .1}
    with pytest.raises(ValueError) as exc:
        run_case(case)
    assert str(exc.value) == "simulation: missing ['end_time'], unsupported []"

    case["simulation"] = {"unused": 1}
    with pytest.raises(ValueError) as exc:
        run_case(case)
    assert str(exc.value) == (
        "simulation: missing ['dt', 'end_time'], unsupported ['unused']"
    )


def test_explicit_ramp_time_remains_distinct_from_omitted_default():
    case = {
        "simulation": {"dt": .1, "end_time": 0},
        "wave": {"type": "none"},
        "bodies": [{"hydro_file": HYDRO}],
        "constraint": {"kind": "heave"},
    }
    assert run_case(case).time.tolist() == [0]

    case["simulation"]["ramp_time"] = 100
    with pytest.raises(ValueError, match="ramp_time is inapplicable to no-wave free decay"):
        run_case(case)

"""Check configurable PTO attachment geometry and public case outputs."""

import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest
from scipy.integrate import trapezoid

from wecsim import NoWave, WEC, WorldPoint
from wecsim.dynamics.ptoConnections import build_linear_ptos
from wecsim.dynamics.caseDynamics import run_case
from wecsim.dynamics.linearCoordinates import build_coordinate_maps


ROOT = Path(__file__).resolve().parents[1]
RM3 = ROOT / "examples/data/rm3.h5"
SPHERE = ROOT / "examples/data/sphere.h5"
EXAMPLE = ROOT / "examples/configurable_rm3_pto.json"


def test_two_body_attachment_offsets_create_pitch_moment():
    first = np.array([
        [0, 0, 0], [0, 0, 0], [1, 0, 0],
        [0, 0, 0], [0, 0, 1], [0, 0, 0],
    ])
    second = np.array([
        [0, 0, 0], [0, 0, 0], [0, 1, 0],
        [0, 0, 0], [0, 0, 0], [0, 0, 0],
    ])
    specs = [{
        "name": "main", "kind": "linear_actuator",
        "from": {"body": "float", "point": [2, 0, 0]},
        "to": {"body": "spar", "point": [-1, 0, 0]},
        "axis": [0, 0, 1], "stiffness": 1000, "damping": 100,
        "equilibrium_position": 0.1,
    }]
    connections, stiffness, damping, bias = build_linear_ptos(
        specs, [first, second], [np.zeros(3), np.zeros(3)],
        ["float", "spar"],
    )
    jacobian = np.array([-1, 1, 2])
    np.testing.assert_allclose(connections[0].stroke_jacobian, jacobian)
    np.testing.assert_allclose(stiffness, 1000 * np.outer(jacobian, jacobian))
    np.testing.assert_allclose(damping, 100 * np.outer(jacobian, jacobian))
    np.testing.assert_allclose(bias, 100 * jacobian)
    outputs = dict(connections[0].outputs(np.zeros((1, 3)), np.zeros((1, 3))))
    np.testing.assert_allclose(outputs["pto_main_force"], [100])
    np.testing.assert_allclose(outputs["pto_main_absorbed_power"], [0])

    # An exact small-pitch attachment displacement gives this Jacobian at
    # the neutral pose. The same Jacobian must map line force back to body
    # coordinates, so virtual work and damper dissipation agree.
    def exact_stroke(q):
        return q[1] - q[0] + 2 * np.sin(q[2])

    step = 1e-6
    finite_difference = np.array([
        (exact_stroke(step * np.eye(3)[axis])
         - exact_stroke(-step * np.eye(3)[axis])) / (2 * step)
        for axis in range(3)
    ])
    np.testing.assert_allclose(finite_difference, jacobian, rtol=0, atol=1e-9)
    position = np.array([0.12, -0.03, 0.04])
    speed = np.array([0.2, -0.3, 0.1])
    line_force = connections[0].force(position[None], speed[None])[0]
    generalized_force = -stiffness @ position - damping @ speed + bias
    np.testing.assert_allclose(generalized_force, jacobian * line_force,
                               rtol=0, atol=1e-12)
    output = dict(connections[0].outputs(position[None], speed[None]))
    stroke = output["pto_main_stroke"][0]
    stroke_speed = output["pto_main_velocity"][0]
    dissipated = output["pto_main_absorbed_power"][0]
    np.testing.assert_allclose(
        generalized_force @ speed + 1000 * (stroke - 0.1) * stroke_speed
        + dissipated, 0, rtol=0, atol=1e-12,
    )


def test_world_anchor_sets_default_axis():
    mapping = np.zeros((6, 1))
    mapping[2, 0] = 1
    specs = [{
        "name": "anchor", "kind": "linear_actuator",
        "from": {"body": 1, "point": [0, 0, 0]},
        "to": {"ground": [0, 0, 2]}, "stiffness": 1000,
    }]
    connections, _, _, _ = build_linear_ptos(
        specs, [mapping], [np.zeros(3)], ["body1"],
    )
    np.testing.assert_allclose(connections[0].stroke_jacobian, [-1])
    np.testing.assert_allclose(
        connections[0].force(np.array([[0.25]]), np.array([[0]])),
        [250],
    )


def test_opposite_body_local_ptos_superpose_as_one_passive_pto():
    """Opposite stroke signs must still give the same generalized force."""
    def run(two_ptos):
        wec = WEC("Sphere PTO superposition")
        body = wec.body("sphere", SPHERE)
        wec.coordinate("heave", body.move("heave"))
        if two_ptos:
            # The upper body-to-ground stroke is -heave; lower ground-to-body
            # stroke is +heave. Their spring biases sum to -2,000 N.
            wec.pto("upper", body.at(0, 0, 0), WorldPoint(0, 0, 5),
                    stiffness=10_000, damping=200_000, pretension=-1_000)
            wec.pto("lower", WorldPoint(0, 0, -5), body.at(0, 0, 0),
                    stiffness=20_000, damping=300_000,
                    equilibrium_position=-.05)
        else:
            wec.pto("equivalent", WorldPoint(0, 0, -5), body.at(0, 0, 0),
                    stiffness=30_000, damping=500_000,
                    equilibrium_position=-1 / 15)
        return wec.run(NoWave(), dt=.01, end_time=2,
                       radiation_memory=.2,
                       initial_coordinate={"heave": .5})

    combined = run(True)
    equivalent = run(False)
    upper, lower = combined.ptos["upper"], combined.ptos["lower"]
    single = equivalent.ptos["equivalent"]
    assert combined.case["ptos"][0]["from"] == {
        "body": "sphere", "point": [0, 0, 0],
    }
    assert combined.case["ptos"][0]["pretension"] == -1_000
    np.testing.assert_allclose(upper.stroke, -lower.stroke, rtol=0, atol=1e-12)
    np.testing.assert_allclose(lower.stroke, single.stroke, rtol=0, atol=1e-12)
    np.testing.assert_allclose(-upper.force + lower.force, single.force,
                               rtol=0, atol=1e-9)
    np.testing.assert_allclose(upper.absorbed_power + lower.absorbed_power,
                               single.absorbed_power, rtol=0, atol=1e-8)
    np.testing.assert_allclose(combined.bodies["sphere"].position,
                               equivalent.bodies["sphere"].position,
                               rtol=0, atol=1e-10)
    np.testing.assert_allclose(combined.bodies["sphere"].velocity,
                               equivalent.bodies["sphere"].velocity,
                               rtol=0, atol=1e-10)
    assert np.ptp(combined.bodies["sphere"].position[:, 2]) > .1
    assert np.max(single.absorbed_power) > 0


def test_public_coordinate_order_does_not_change_coupled_body_or_pto_motion():
    """Named initial states and body-local PTO geometry are order independent."""
    def run(order):
        wec = WEC("Sphere coordinate order")
        body = wec.body("sphere", SPHERE)
        for name in order:
            wec.coordinate(name, body.move(name))
        wec.pto(
            "heave_damper", WorldPoint(0, 0, -5), body.at(0, 0, 0),
            stiffness=25_000, damping=50_000,
        )
        return wec.run(
            NoWave(), dt=.02, end_time=1, radiation_memory=.2,
            initial_coordinate={"surge": .1, "heave": .3},
            initial_speed={"surge": -.02, "heave": .05},
        )

    surge_first = run(("surge", "heave"))
    heave_first = run(("heave", "surge"))
    for field in ("position", "velocity"):
        np.testing.assert_array_equal(
            getattr(surge_first.bodies["sphere"], field),
            getattr(heave_first.bodies["sphere"], field),
        )
    for field in ("stroke", "velocity", "force", "absorbed_power"):
        np.testing.assert_array_equal(
            getattr(surge_first.ptos["heave_damper"], field),
            getattr(heave_first.ptos["heave_damper"], field),
        )
    assert np.ptp(surge_first.bodies["sphere"].position[:, 0]) > .01
    assert np.ptp(surge_first.bodies["sphere"].position[:, 2]) > .01
    assert np.max(surge_first.ptos["heave_damper"].absorbed_power) > 0


def test_position_feedback_matches_spring_motion_but_accounts_for_its_work():
    """Equivalent forces need not have the same reported absorbed energy."""
    spring = 25_000.0
    damper = 50_000.0

    def run(*, feedback):
        wec = WEC("Sphere heave with configurable PTO")
        body = wec.body("sphere", SPHERE)
        wec.coordinate("heave", body.move("heave"))
        wec.pto(
            "main", WorldPoint(0, 0, -5), body.at(0, 0, 0),
            damping=damper,
            position_control_gain=spring if feedback else 0,
            stiffness=0 if feedback else spring,
        )
        return wec.run(
            NoWave(), dt=.01, end_time=2, radiation_memory=.2,
            initial_coordinate={"heave": .5},
        )

    passive = run(feedback=False)
    active = run(feedback=True)
    passive_pto = passive.ptos["main"]
    active_pto = active.ptos["main"]
    for field in ("position", "velocity"):
        np.testing.assert_allclose(
            getattr(passive.bodies["sphere"], field),
            getattr(active.bodies["sphere"], field), rtol=0, atol=1e-10,
        )
    for field in ("stroke", "velocity", "force"):
        np.testing.assert_allclose(
            getattr(passive_pto, field), getattr(active_pto, field),
            rtol=0, atol=1e-10,
        )
    assert np.ptp(active_pto.stroke) > .1
    np.testing.assert_allclose(
        passive_pto.absorbed_power, damper * passive_pto.velocity**2,
        rtol=0, atol=1e-10,
    )
    np.testing.assert_allclose(
        active_pto.control_force, -spring * active_pto.stroke,
        rtol=0, atol=1e-10,
    )
    np.testing.assert_allclose(
        active_pto.absorbed_power - passive_pto.absorbed_power,
        active_pto.control_absorbed_power, rtol=0, atol=1e-9,
    )
    # The reactive controller exchanges the energy a passive spring stores.
    controller_work = trapezoid(active_pto.control_absorbed_power,
                                x=active.time)
    spring_energy_change = spring * (
        active_pto.stroke[-1]**2 - active_pto.stroke[0]**2
    ) / 2
    assert controller_work < 0
    np.testing.assert_allclose(controller_work, spring_energy_change,
                               rtol=1e-4, atol=0)


def test_named_pitch_pivot_moves_body_center():
    constraint = {"coordinates": [{
        "name": "pitch", "motions": [{
            "body": "float", "dof": "pitch",
            "pivot": {"world": [0, 0, 0]},
        }],
    }]}
    maps, names = build_coordinate_maps(
        constraint, [{}], ["float"], [np.array([0, 0, -0.72])],
    )
    assert names == ("pitch",)
    np.testing.assert_allclose(maps[0][:, 0], [-0.72, 0, 0, 0, 1, 0])
    constraint["coordinates"][0]["motions"][0]["pivot"] = {
        "point": [0, 0, 0.72],
    }
    local_maps, _ = build_coordinate_maps(
        constraint, [{}], ["float"], [np.array([0, 0, -0.72])],
    )
    np.testing.assert_allclose(local_maps[0], maps[0])


def test_connection_case_writes_stroke_force_and_power(tmp_path):
    case = {
        "simulation": {"dt": 0.1, "end_time": 0.1, "ramp_time": 0},
        "wave": {"type": "regular", "height": 0, "period": 8},
        "bodies": [
            {"name": "float", "hydro_file": str(RM3),
             "inertia": [0, 1e6, 0]},
            {"name": "spar", "hydro_file": str(RM3)},
        ],
        "constraint": {
            "kind": "linear_subspace",
            "coordinates": [
                {"name": "float_heave", "motions": [
                    {"body": "float", "dof": "heave"}]},
                {"name": "spar_heave", "motions": [
                    {"body": "spar", "dof": "heave"}]},
                {"name": "float_pitch", "motions": [
                    {"body": "float", "dof": "pitch",
                     "pivot": {"world": [0, 0, 0]}}]},
            ],
            "initial_coordinate": {"float_pitch": 0},
        },
        "ptos": [{
            "name": "main", "kind": "linear_actuator",
            "from": {"body": "float", "point": [2, 0, 0.72]},
            "to": {"body": "spar", "point": [-1, 0, 21.29]},
            "axis": [0, 0, 1], "stiffness": 1000,
            "equilibrium_position": 0.1,
        }],
    }
    case_file = tmp_path / "case.json"
    output = tmp_path / "motion.csv"
    case_file.write_text(json.dumps(case), encoding="utf-8")
    result = subprocess.run(
        [sys.executable, "-m", "wecsim",
         str(case_file), "--output", str(output)],
        cwd=ROOT, capture_output=True, text=True, check=False,
    )
    assert result.returncode == 0, result.stderr
    values = np.loadtxt(output, delimiter=",", skiprows=1)
    columns = output.read_text(encoding="utf-8").splitlines()[0].split(",")
    data = dict(zip(columns, values.T))
    np.testing.assert_allclose(data["pto_main_stroke"][0], 0, atol=1e-12)
    np.testing.assert_allclose(data["pto_main_force"][0], 100, atol=1e-12)
    np.testing.assert_allclose(data["pto_coordinate1_force"][0], -100, atol=1e-12)
    np.testing.assert_allclose(data["pto_coordinate2_force"][0], 100, atol=1e-12)
    np.testing.assert_allclose(data["pto_coordinate3_force"][0], 200, atol=1e-12)
    assert "coordinate_float_pitch_position" in data
    assert "pto_main_absorbed_power" in data
    assert np.isfinite(values).all()


def test_connection_rejects_unknown_body():
    mapping = np.zeros((6, 1))
    mapping[2, 0] = 1
    specs = [{
        "name": "bad", "kind": "linear_actuator",
        "from": {"body": "missing", "point": [0, 0, 0]},
        "to": {"ground": [0, 0, 2]}, "stiffness": 1000,
    }]
    with pytest.raises(ValueError, match="identify a body"):
        build_linear_ptos(specs, [mapping], [np.zeros(3)], ["body1"])


@pytest.mark.parametrize("change, message", [
    ({"to": {"ground": [0, 0, 0]}}, "noncoincident endpoints"),
    ({"to": {"ground": [2, 0, 0]}}, "no motion in its force axis"),
    ({"from": {"body": True, "point": [0, 0, 0]}},
     "identify a body"),
    ({"from": {"body": 1, "point": [0, 0, float("nan")]}},
     "finite three-vector"),
    ({"damping": -1}, "outside its supported range"),
    ({"pretension": 100}, "pretension needs positive stiffness"),
    ({"equilibrium_position": 0.1}, "equilibrium needs positive stiffness"),
], ids=["coincident-points", "uncoupled-axis", "boolean-body",
        "nonfinite-attachment", "negative-damping", "unsupported-pretension",
        "unsupported-equilibrium"])
def test_pto_configuration_rejects_nonphysical_force_paths(change, message):
    """A user-defined actuator must have a finite, movable line of action."""
    mapping = np.zeros((6, 1))
    mapping[2, 0] = 1  # The body can heave, but cannot surge.
    spec = {
        "name": "main", "kind": "linear_actuator",
        "from": {"body": 1, "point": [0, 0, 0]},
        "to": {"ground": [0, 0, 2]},
        "damping": 100,
        **change,
    }
    with pytest.raises(ValueError, match=message):
        build_linear_ptos([spec], [mapping], [np.zeros(3)], ["body1"])


def test_pto_configuration_rejects_duplicate_output_names():
    mapping = np.zeros((6, 1))
    mapping[2, 0] = 1
    spec = {
        "name": "main", "kind": "linear_actuator",
        "from": {"body": 1, "point": [0, 0, 0]},
        "to": {"ground": [0, 0, 2]},
        "damping": 100,
    }
    with pytest.raises(ValueError, match="unique identifiers"):
        build_linear_ptos([spec, spec], [mapping], [np.zeros(3)], ["body1"])


def test_configurable_example_runs_with_named_coordinates():
    case = json.loads(EXAMPLE.read_text(encoding="utf-8"))
    case["simulation"]["end_time"] = 0.2
    response = run_case(case, base_dir=EXAMPLE.parent)
    outputs = dict(response.extra_outputs)
    assert response.body_position.shape == (3, 2, 6)
    assert "coordinate_shared_pitch_position" in outputs
    assert "pto_main_stroke" in outputs
    assert np.isfinite(outputs["pto_main_force"]).all()

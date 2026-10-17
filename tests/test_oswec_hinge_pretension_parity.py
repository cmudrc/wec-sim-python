"""Pair a configured OSWEC hinge spring and pretension with pinned WEC-Sim."""

import os
from pathlib import Path

import numpy as np

from wecsim import RegularWave, WEC, WorldPoint


APPLICATIONS = os.environ.get("WEC_SIM_APPLICATIONS_DIR")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests/fixtures/matlab/oswec_hinge_pretension/source.npz"
PINNED_HYDRO = ROOT / "tests/fixtures/matlab/oswec_multi_wave/oswec.h5"


def _source():
    """Prefer fresh MATLAB CSVs; keep the same full-duration gate in CI."""
    if APPLICATIONS and REFERENCE:
        folder = Path(REFERENCE)
        if (folder / "OSWEC_Nonhydro_Pretension_parameters.csv").is_file():
            paths = {
                "flap": "OSWEC_Nonhydro_Pretension_body1.csv",
                "base": "OSWEC_Nonhydro_Pretension_body2.csv",
                "pto": "OSWEC_Nonhydro_Pretension_pto1.csv",
                "published_flap": "OSWEC_Nonhydro_Nonhydro_body1.csv",
                "parameters": "OSWEC_Nonhydro_Pretension_parameters.csv",
            }
            records = {name: np.loadtxt(folder / filename, delimiter=",")
                       for name, filename in paths.items()}
            hydro = (Path(APPLICATIONS)
                     / "_Common_Input_Files/OSWEC/hydroData/oswec.h5")
            return records, hydro.resolve()
    with np.load(FIXTURE, allow_pickle=False) as saved:
        records = {name: saved[name].copy() for name in saved.files}
    return records, PINNED_HYDRO


def _max_error(actual, expected, limit, label):
    assert actual.shape == expected.shape, label
    assert np.isfinite(actual).all() and np.isfinite(expected).all(), label
    error = float(np.max(np.abs(actual - expected)))
    print(f"{label}: {error:.6g} (gate {limit:.6g})")
    assert error < limit, label


def test_configured_hinge_spring_and_pretension_against_matlab():
    source, hydro = _source()
    parameters = source["parameters"]
    np.testing.assert_allclose(
        parameters, [50_000, 12_000, -5_000, 0.1], rtol=0, atol=1e-12,
    )
    bodies = [source["flap"], source["base"]]
    published_flap = source["published_flap"]
    pto = source["pto"]
    assert all(values.shape == (4001, 25) for values in (*bodies, pto))
    assert published_flap.shape == (4001, 25)
    # The pinned rotational mask treats equilibriumPosition as degrees and
    # applies its spring offset opposite to the reported positive pitch.
    # This source pretension therefore produces -87.266 N m at zero pitch,
    # not the +5,000 N m implied by a physical N m pretension input.
    source_equilibrium_rad = np.deg2rad(parameters[3])
    equivalent_pretension = -parameters[2] * np.pi / 180
    np.testing.assert_allclose(
        equivalent_pretension,
        parameters[0] * source_equilibrium_rad,
        rtol=0, atol=1e-12,
    )
    np.testing.assert_allclose(
        pto[0, 17], parameters[2] * np.pi / 180,
        rtol=0, atol=1e-10,
    )
    assert abs(pto[0, 17] + parameters[2]) > 1_000

    wec = WEC("OSWEC configured hinge PTO")
    flap = wec.body("flap", hydro, mass=127_000,
                    inertia=(0, 1.85e6, 0))
    base = wec.fixed_body(
        "base", center_gravity=(0, 0, -10.9), mass=999,
        inertia=(1, 1, 1),
    )
    wec.fixed_hinge(
        flap, base, location=WorldPoint(0, 0, -10),
        pto_location=WorldPoint(0, 0, -8.9), damping=12_000,
        stiffness=50_000, pretension=equivalent_pretension,
    )
    result = wec.run(
        RegularWave(2.5, 8), dt=0.1, end_time=400, ramp_time=100,
    )
    np.testing.assert_allclose(
        result.case["pto"]["equilibrium_position"],
        -source_equilibrium_rad, rtol=0, atol=1e-15,
    )
    _max_error(result.time, pto[:, 0], 1e-10, "time")
    for name, expected in zip(("flap", "base"), bodies):
        _max_error(result.bodies[name].position,
                   expected[:, 1:7],
                   0.03 if name == "flap" else 1e-10,
                   f"{name} position")
        _max_error(result.bodies[name].velocity,
                   expected[:, 7:13],
                   0.02 if name == "flap" else 1e-10,
                   f"{name} velocity")

    history = result.ptos["hinge"]
    _max_error(history.stroke, pto[:, 5], 0.006, "PTO angle")
    _max_error(history.velocity, pto[:, 11], 0.005, "PTO angular speed")
    _max_error(history.force, pto[:, 17], 350, "PTO torque")
    _max_error(history.force,
               -12_000 * history.velocity
               - 50_000 * (history.stroke + source_equilibrium_rad),
               1e-8, "Python PTO torque law")
    _max_error(pto[:, 17],
               -12_000 * pto[:, 11]
               - 50_000 * (pto[:, 5] + source_equilibrium_rad),
               1e-5, "MATLAB PTO torque law")
    _max_error(pto[:, 23], pto[:, 17] * pto[:, 11],
               1e-5, "MATLAB PTO mechanical power")
    _max_error(-history.force * history.velocity, -pto[:, 23],
               120, "PTO mechanical power")

    # Paired absolute motion errors include the original fixed-hinge solver
    # difference. Compare each implementation's spring-induced motion to its
    # own no-PTO control to isolate the response to this configured spring.
    control = WEC("OSWEC published no-PTO control")
    control_flap = control.body("flap", hydro, mass=127_000,
                                inertia=(0, 1.85e6, 0))
    control_base = control.fixed_body(
        "base", center_gravity=(0, 0, -10.9), mass=999,
        inertia=(1, 1, 1),
    )
    control.fixed_hinge(
        control_flap, control_base, location=WorldPoint(0, 0, -10),
        pto_location=WorldPoint(0, 0, -8.9),
    )
    control_result = control.run(
        RegularWave(2.5, 8), dt=0.1, end_time=400, ramp_time=100,
    )
    for dof, label, position_limit, speed_limit, min_effect in (
        (0, "surge", 0.00015, 0.0001, 0.002),
        (2, "heave", 0.00008, 0.00007, 0.0007),
        (4, "pitch", 0.00003, 0.00003, 0.0005),
    ):
        source_effect = bodies[0][:, 1 + dof] - published_flap[:, 1 + dof]
        python_effect = (result.bodies["flap"].position[:, dof]
                         - control_result.bodies["flap"].position[:, dof])
        assert np.max(np.abs(source_effect)) > min_effect
        _max_error(python_effect, source_effect, position_limit,
                   f"spring-induced {label} position")
        source_speed_effect = (bodies[0][:, 7 + dof]
                               - published_flap[:, 7 + dof])
        python_speed_effect = (result.bodies["flap"].velocity[:, dof]
                               - control_result.bodies["flap"].velocity[:, dof])
        _max_error(python_speed_effect, source_speed_effect, speed_limit,
                   f"spring-induced {label} speed")


def test_general_coordinate_hinge_and_rotational_pto_against_matlab():
    """The public coordinate route must use the same nonlinear hinge motion."""
    source, hydro = _source()
    flap_source = source["flap"]
    pto_source = source["pto"]
    wec = WEC("OSWEC coordinate hinge")
    flap = wec.body("flap", hydro, mass=127_000,
                    inertia=(0, 1.85e6, 0))
    pitch = wec.coordinate(
        "pitch", flap.move("pitch", pivot=WorldPoint(0, 0, -8.9)),
    )
    wec.rotational_pto(
        "hinge", pitch, damping=12_000, stiffness=50_000,
        pretension=5_000 * np.pi / 180,
    )
    result = wec.run(
        RegularWave(2.5, 8), dt=0.1, end_time=400, ramp_time=100,
    )
    assert flap_source.shape == pto_source.shape == (4001, 25)
    np.testing.assert_allclose(
        result.case["ptos"][0]["equilibrium_position"],
        -np.deg2rad(0.1), rtol=0, atol=1e-15,
    )
    _max_error(result.bodies["flap"].position, flap_source[:, 1:7],
               0.03, "coordinate-hinge flap position")
    _max_error(result.bodies["flap"].velocity, flap_source[:, 7:13],
               0.02, "coordinate-hinge flap velocity")
    history = result.ptos["hinge"]
    _max_error(history.stroke, pto_source[:, 5],
               0.006, "coordinate-hinge PTO angle")
    _max_error(history.velocity, pto_source[:, 11],
               0.005, "coordinate-hinge PTO speed")
    _max_error(history.force, pto_source[:, 17],
               350, "coordinate-hinge PTO torque")
    _max_error(-history.force * history.velocity, -pto_source[:, 23],
               120, "coordinate-hinge PTO mechanical power")
    lever = result.bodies["flap"].position[0, 2] + 8.9
    angle = result.coordinates["pitch"].position
    _max_error(result.bodies["flap"].position[:, 0],
               lever * np.sin(angle), 1e-10,
               "exact pivot surge position")
    _max_error(result.bodies["flap"].position[:, 2],
               -8.9 + lever * np.cos(angle), 1e-10,
               "exact pivot heave position")

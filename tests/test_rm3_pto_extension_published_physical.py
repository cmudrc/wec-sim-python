"""Portable physical gates for the unchanged RM3 PTO Extension inputs."""

import hashlib
import os
from pathlib import Path

import h5py
import numpy as np
import pytest

from wecsim.dynamics.caseDynamics import run_case


ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests/fixtures/matlab/rm3_pto_extension_published_physical"
SOURCE = FIXTURE / "source.npz"
HYDRO = FIXTURE / "rm3.h5"
STEP = ROOT / "tests/fixtures/matlab/rm3_pto_extension_step.npz"
FRESH_SOURCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
FRESH_APPLICATIONS = os.environ.get("WEC_SIM_APPLICATIONS_DIR")

# Retain the limits declared in test_rm3_pto_extension_parity.py before the
# source-step investigation. Only the qualifying solver changes here: this
# test uses the public physical implicit-mass path, not source-style ode4.
CASES = (
    ("float", {"float_heave": 5}, 0, 0.090, 0.210, 0.050),
    ("spar", {"spar_heave": -5}, 1, 0.018, 0.008, 0.008),
)


def _bound(actual, expected, limit, label):
    assert actual.shape == expected.shape, label
    assert np.isfinite(actual).all() and np.isfinite(expected).all(), label
    error = float(np.max(np.abs(actual - expected)))
    print(f"{label}: {error:.9g} (limit {limit:.9g})")
    assert error < limit, label
    return error


def _case(coordinate, *, numerical=False):
    simulation = {
        "dt": 0.1, "end_time": 30, "radiation_memory": 60,
        "added_mass_scheme": "simulink_delay" if numerical else "implicit",
    }
    if numerical:
        simulation["solver"] = "wecsim_ode4"
    return {
        "simulation": simulation,
        "wave": {"type": "none"},
        "bodies": [
            {"hydro_file": str(HYDRO), "hydro_body": 1,
             "mass": "equilibrium", "pitch_inertia": 21_306_090.66},
            {"hydro_file": str(HYDRO), "hydro_body": 2,
             "mass": "equilibrium", "pitch_inertia": 94_407_091.24},
        ],
        "constraint": {"kind": "floating_joint", "location": [0, 0, 0],
                       "initial_coordinate": coordinate},
        "pto": {"kind": "relative_heave", "damping": 0},
    }


def test_published_physical_fixture_matches_pinned_source():
    assert hashlib.sha256(HYDRO.read_bytes()).hexdigest() == (
        "dac35aa42321e1ade04e8935b38680d7b70e7b247edac9b5b800d01142fcc259"
    )
    assert hashlib.sha256(SOURCE.read_bytes()).hexdigest() == (
        "535eaabbf0ec68949ffffff215548003696a0ad99029214a396b37d6fa67e00c"
    )
    active = (
        "properties/cg", "properties/disp_vol",
        "hydro_coeffs/linear_restoring_stiffness",
        "hydro_coeffs/added_mass/inf_freq",
        "hydro_coeffs/radiation_damping/impulse_response_fun/K",
        "hydro_coeffs/radiation_damping/impulse_response_fun/t",
    )
    if FRESH_APPLICATIONS:
        apps = Path(FRESH_APPLICATIONS)
        for name in ("float", "spar"):
            input_file = apps / "RM3_PTO_Extension" / name / "wecSimInputFile.m"
            if input_file.is_file():
                contents = input_file.read_text()
                for setting in (
                    "simu.endTime = 30;", "simu.solver = 'ode4';",
                    "simu.dt = 0.1;", "waves = waveClass('noWaveCIC');",
                    "pto(1).stiffness = 0;", "pto(1).damping = 0;",
                    "pto(1).extension.PositionTargetValue = 5;",
                ):
                    assert setting in contents, (name, setting)
    with h5py.File(HYDRO) as source, h5py.File(
        ROOT / "examples/data/rm3_end_stops.h5"
    ) as bundled:
        for body in ("body1", "body2"):
            for suffix in active:
                key = f"{body}/{suffix}"
                np.testing.assert_array_equal(source[key], bundled[key])
        if FRESH_APPLICATIONS:
            fresh_hydro = (Path(FRESH_APPLICATIONS) /
                           "_Common_Input_Files/RM3/hydroData/rm3.h5")
            with h5py.File(fresh_hydro) as fresh:
                for body in ("body1", "body2"):
                    for suffix in active:
                        key = f"{body}/{suffix}"
                        if suffix.endswith("impulse_response_fun/K"):
                            np.testing.assert_allclose(
                                source[key], fresh[key], rtol=0, atol=1e-10,
                            )
                        else:
                            np.testing.assert_array_equal(source[key], fresh[key])

    with np.load(SOURCE, allow_pickle=False) as saved, \
            np.load(STEP, allow_pickle=False) as prior:
        assert set(saved.files) == {
            f"{name}_{part}" for name in ("float", "spar")
            for part in ("body1", "body2", "pto1")
        }
        for name in ("float", "spar"):
            body1, body2, pto = (
                saved[f"{name}_{part}"] for part in ("body1", "body2", "pto1")
            )
            for array in (body1, body2, pto):
                assert array.shape == (301, 25) and np.isfinite(array).all()
                np.testing.assert_allclose(
                    array[:, 0], np.arange(301) * .1, rtol=0, atol=1e-10,
                )
            assert np.max(np.abs(pto[:, 15])) < 1e-7
            assert np.max(np.abs(pto[:, 21])) < 1e-7
            compact = np.column_stack((
                body1[:, 0], body1[:, 3], body1[:, 9],
                body2[:, 3], body2[:, 9],
                pto[:, 3], pto[:, 9], pto[:, 15], pto[:, 21],
            ))
            np.testing.assert_array_equal(compact, prior[f"{name}_0p1"])
            if FRESH_SOURCE:
                for part, archived in (("body1", body1),
                                       ("body2", body2), ("pto1", pto)):
                    fresh = np.loadtxt(
                        Path(FRESH_SOURCE) /
                        f"RM3_PTO_Extension_{name}_{part}.csv",
                        delimiter=",", ndmin=2,
                    )
                    np.testing.assert_allclose(fresh, archived, rtol=1e-10,
                                               atol=1e-6)


@pytest.mark.parametrize(
    "name,coordinate,active_body,heave_limit,speed_limit,control_gap",
    CASES,
)
def test_unchanged_published_input_with_independent_physical_motion(
    name, coordinate, active_body, heave_limit, speed_limit, control_gap,
):
    # MATLAB states, accelerations, forces, and solver steps are read only
    # after the independent physical run has completed.
    response = run_case(_case(coordinate))
    assert response.wave_elevation is None
    assert response.body_position.shape == response.body_velocity.shape == (
        301, 2, 6,
    )
    with np.load(SOURCE, allow_pickle=False) as saved:
        source = {part: saved[f"{name}_{part}"].copy()
                  for part in ("body1", "body2", "pto1")}

    for index in range(2):
        expected = source[f"body{index + 1}"]
        _bound(response.time, expected[:, 0], 1e-10,
               f"{name} body{index + 1} time")
        for axis in range(6):
            if axis == 2:
                limits = ((heave_limit, speed_limit) if index == active_body
                          else (1e-8, 1e-8))
            elif axis == 0:
                limits = (3e-5, 1e-5)
            elif axis == 4:
                limits = (1e-6, 1e-6)
            else:
                limits = (1e-9, 1e-9)
            _bound(response.body_position[:, index, axis],
                   expected[:, 1 + axis], limits[0],
                   f"{name} body{index + 1} axis {axis} position")
            _bound(response.body_velocity[:, index, axis],
                   expected[:, 7 + axis], limits[1],
                   f"{name} body{index + 1} axis {axis} speed")

    with h5py.File(HYDRO) as hydro:
        center_gap = (
            np.asarray(hydro["body1/properties/cg"]).ravel()[2]
            - np.asarray(hydro["body2/properties/cg"]).ravel()[2]
        )
    pitch = response.body_position[:, 0, 4]
    pitch_speed = response.body_velocity[:, 0, 4]
    cosine = np.cos(pitch)
    stroke = ((response.body_position[:, 0, 2]
               - response.body_position[:, 1, 2]) / cosine - center_gap)
    speed = ((response.body_velocity[:, 0, 2]
              - response.body_velocity[:, 1, 2]
              + (center_gap + stroke) * np.sin(pitch) * pitch_speed) / cosine)
    pto = source["pto1"]
    _bound(response.time, pto[:, 0], 1e-10, f"{name} PTO time")
    assert abs(pto[0, 3] - 5) < 1e-10
    assert abs(stroke[0] - 5) < 1e-8
    assert abs(speed[0]) < 1e-8
    _bound(stroke, pto[:, 3], heave_limit, f"{name} PTO stroke")
    _bound(speed, pto[:, 9], speed_limit, f"{name} PTO speed")
    _bound(response.pto_force, pto[:, 15], 1e-7, f"{name} PTO zero force")
    _bound(response.pto_force * speed, pto[:, 21], 1e-7,
           f"{name} PTO zero power")

    # Strict negative control: the opt-in source-numerical ode4 path nearly
    # reproduces the archived trace, while the qualifying implicit trajectory
    # is distinct. A switch of the qualifying run to that route fails here.
    numerical = run_case(_case(coordinate, numerical=True))
    reference_heave = source[f"body{active_body + 1}"][:, 3]
    _bound(numerical.body_position[:, active_body, 2], reference_heave,
           1e-8, f"{name} source-numerical control")
    separation = float(np.max(np.abs(
        response.body_position[:, active_body, 2]
        - numerical.body_position[:, active_body, 2]
    )))
    print(f"{name} physical/source-numerical separation: {separation:.9g}")
    assert separation > control_gap

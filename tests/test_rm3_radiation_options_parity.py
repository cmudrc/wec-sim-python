"""Compare published RM3 radiation options with paired MATLAB trajectories."""

import os
from pathlib import Path

import h5py
import numpy as np
import pytest
from scipy.interpolate import CubicSpline
from scipy.signal import fftconvolve

from wecsim.dynamics.caseDynamics import run_case
from wecsim.hydrodynamics.radiation import replay_rm3_fitted_radiation


APPLICATIONS = os.environ.get("WEC_SIM_APPLICATIONS_DIR")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
FIXTURE = Path(__file__).parent / "fixtures/matlab/rm3_radiation_options"
if bool(APPLICATIONS) != bool(REFERENCE):
    raise ValueError("provide both RM3 radiation source paths or neither")
HYDRO = (Path(APPLICATIONS) / "_Common_Input_Files/RM3/hydroData/rm3.h5"
         if APPLICATIONS else FIXTURE / "rm3.h5")


def _source(mode, part):
    if REFERENCE:
        return np.loadtxt(
            Path(REFERENCE) / f"RM3_Radiation_Options_{mode}_{part}.csv",
            delimiter=",", ndmin=2,
        )
    with np.load(FIXTURE / "source.npz") as source:
        return source[f"{mode}_{part}"]


def _max_error(actual, expected, limit, label):
    assert actual.shape == expected.shape, label
    assert np.isfinite(expected).all(), label
    error = np.max(np.abs(actual - expected))
    assert error < limit, f"{label}: max error {error:.6g} exceeds {limit}"


def _fir_radiation_force(h5, body, velocity, dt=0.1, memory=60):
    """Reconstruct each body's FIR output from source HDF5 and body velocity."""
    path = f"body{body}/hydro_coeffs/radiation_damping/impulse_response_fun"
    kernel = np.asarray(h5[f"{path}/K"])
    kernel_time = np.asarray(h5[f"{path}/t"]).ravel()
    rho = float(np.asarray(h5["simulation_parameters/rho"]).item())
    taps = rho * CubicSpline(kernel_time, kernel, axis=2)(
        np.arange(round(memory / dt) + 1) * dt,
    )
    force = np.zeros_like(velocity)
    for output in range(6):
        for input_dof in range(6):
            force[:, output] += dt * fftconvolve(
                velocity[:, input_dof],
                taps[output, 6 * (body - 1) + input_dof],
                mode="full",
            )[:len(velocity)]
    return force


@pytest.mark.parametrize("mode", ["constant", "convolution", "FIR",
                                  "state_space"])
def test_published_rm3_radiation_options(mode):
    hydro = HYDRO.resolve()
    simulation = {"dt": 0.1, "end_time": 500, "ramp_time": 100}
    if mode != "constant":
        simulation["radiation_memory"] = 60
    if mode == "convolution":
        simulation["added_mass_scheme"] = "simulink_delay"
    if mode == "FIR":
        simulation["radiation_method"] = "fir"
    if mode == "state_space":
        # Reproduce the published nonpassive fit explicitly; convolution
        # remains the default radiation model.
        simulation["radiation_method"] = "state_space"
    case = {
        "simulation": simulation,
        "wave": {"type": "regularCIC" if mode != "constant" else "regular",
                 "height": 2.5, "period": 12},
        "bodies": [
            {"hydro_file": str(hydro), "hydro_body": 1,
             "mass": "equilibrium", "pitch_inertia": 21_306_090.66},
            {"hydro_file": str(hydro), "hydro_body": 2,
             "mass": "equilibrium", "pitch_inertia": 94_407_091.24},
        ],
        "constraint": {"kind": "floating_joint", "location": [0, 0, 0]},
        "pto": {"kind": "relative_heave", "damping": 1_200_000},
    }
    response = run_case(case)
    assert response.body_position.shape == (5001, 2, 6)
    for index in range(2):
        expected = _source(mode, f"body{index + 1}")
        assert expected.shape == (5001, 37)  # includes radiation and added-mass forces
        np.testing.assert_allclose(response.time, expected[:, 0],
                                   rtol=0, atol=1e-10)
        for dof, axis, position_limit, velocity_limit in (
            (0, "surge", 0.035, 0.020),
            (2, "heave", 0.012, 0.0025),
            (4, "pitch", 0.0018, 0.0011),
        ):
            _max_error(response.body_position[:, index, dof],
                       expected[:, 1 + dof], position_limit,
                       f"{mode} body{index + 1} {axis} position")
            _max_error(response.body_velocity[:, index, dof],
                       expected[:, 7 + dof], velocity_limit,
                       f"{mode} body{index + 1} {axis} velocity")
        for dof in (1, 3, 5):
            _max_error(response.body_position[:, index, dof],
                       expected[:, 1 + dof], 1e-9,
                       f"{mode} body{index + 1} stationary position {dof}")
            _max_error(response.body_velocity[:, index, dof],
                       expected[:, 7 + dof], 1e-9,
                       f"{mode} body{index + 1} stationary velocity {dof}")

    with h5py.File(hydro) as h5:
        center_gap = (
            np.asarray(h5["body1/properties/cg"]).ravel()[2]
            - np.asarray(h5["body2/properties/cg"]).ravel()[2]
        )
    pitch = response.body_position[:, 0, 4]
    pitch_speed = response.body_velocity[:, 0, 4]
    cosine = np.cos(pitch)
    stroke = ((response.body_position[:, 0, 2]
               - response.body_position[:, 1, 2]) / cosine - center_gap)
    speed = ((response.body_velocity[:, 0, 2]
              - response.body_velocity[:, 1, 2]
              + (center_gap + stroke) * np.sin(pitch) * pitch_speed) / cosine)
    pto = _source(mode, "pto1")
    np.testing.assert_allclose(response.time, pto[:, 0], rtol=0, atol=1e-10)
    _max_error(stroke, pto[:, 3], 0.0125, f"{mode} PTO stroke")
    _max_error(speed, pto[:, 9], 0.003, f"{mode} PTO speed")
    _max_error(response.pto_force, pto[:, 15], 3500,
               f"{mode} PTO force")
    _max_error(response.pto_force * speed, pto[:, 21], 3000,
               f"{mode} PTO power")

    # A samplewise power limit alone does not bound accumulated production.
    # Check the ramp and two developed-wave intervals as well as the 500 s
    # total, keeping each radiation implementation independently accountable.
    python_power = response.pto_force * speed
    source_power = pto[:, 21]
    for start, stop in ((0, 500), (0, 100), (100, 300), (300, 500)):
        window = (response.time >= start) & (response.time <= stop)
        source_work = np.trapezoid(source_power[window], pto[window, 0])
        python_work = np.trapezoid(python_power[window], response.time[window])
        assert source_work < -1e6
        relative_gap = abs(python_work - source_work) / abs(source_work)
        assert relative_gap < .01, (
            f"{mode} PTO work {start}–{stop} s: "
            f"{relative_gap:.4%} exceeds 1%"
        )

    if mode == "FIR":
        with h5py.File(hydro) as h5:
            for body in (1, 2):
                expected = _source("FIR", f"body{body}")
                calculated = _fir_radiation_force(
                    h5, body, response.body_velocity[:, body - 1, :],
                )
                for dof, name, limit in (
                    (0, "surge", 2_000 if body == 1 else 500),
                    (2, "heave", 1_300 if body == 1 else 150),
                    (4, "pitch", 10_000 if body == 1 else 3_000),
                ):
                    _max_error(calculated[:, dof], expected[:, 25 + dof],
                               limit, f"FIR body{body} {name} radiation force")


def test_matlab_radiation_modes_are_distinct():
    loads = {
        mode: _source(mode, "body1")
        for mode in ("constant", "FIR", "state_space", "convolution")
    }
    for response in loads.values():
        assert response.shape == (5001, 37)
        assert np.isfinite(response).all()
    reference_surge = loads["convolution"][:, 1]
    assert np.max(np.abs(loads["FIR"][:, 1] - reference_surge)) > 0.2
    assert np.max(np.abs(loads["state_space"][:, 1] - reference_surge)) > 1.0


def test_published_state_space_radiation_force_on_source_motion():
    """Check the source fitted force law without making it the default solver."""
    hydro = HYDRO.resolve()
    traces = [_source("state_space", f"body{body}") for body in (1, 2)]
    assert all(trace.shape == (5001, 37) for trace in traces)
    time = traces[0][:, 0]
    np.testing.assert_allclose(traces[1][:, 0], time, rtol=0, atol=1e-10)
    velocity = np.stack([trace[:, 7:13] for trace in traces], axis=1)
    expected = np.stack([trace[:, 25:31] for trace in traces], axis=1)
    calculated = replay_rm3_fitted_radiation(hydro, time, velocity)
    error = np.max(np.abs(calculated - expected), axis=0)
    peak = np.max(np.abs(expected), axis=0)
    assert np.all(error < 3e-4 * np.maximum(peak, 3)), (error, peak)

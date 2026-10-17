"""Pair the published Sphere Mean_Drift application with pinned MATLAB."""

import os
import json
from pathlib import Path

import h5py
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.signal import fftconvolve

from wecsim import RegularCICWave, WEC
from wecsim.interface.cli import main as cli_main


APPLICATIONS = os.environ.get("WEC_SIM_APPLICATIONS_DIR")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
PORTABLE = Path(__file__).parent / "fixtures/matlab/sphere_mean_drift"


def _paired_data():
    """Use fresh source output in MATLAB CI and the same archived case elsewhere."""
    if APPLICATIONS or REFERENCE:
        if not (APPLICATIONS and REFERENCE):
            raise ValueError("both Applications and MATLAB output are required")
        root = Path(REFERENCE)
        # The CLI test writes its case to tmp_path. Keep the source checkout
        # address stable when CI supplies a relative Applications directory.
        hydro = Path(APPLICATIONS).resolve() / "Mean_Drift/hydroData/sphere.h5"
        body = np.loadtxt(root / "SPHERE_MEAN_DRIFT_Mean_Drift_body1.csv",
                          delimiter=",")
        wave = np.loadtxt(root / "SPHERE_MEAN_DRIFT_wave.csv", delimiter=",")
        coefficients = np.loadtxt(
            root / "SPHERE_MEAN_DRIFT_excitation_coefficients.csv",
            delimiter=",",
        )
        mass = np.loadtxt(root / "SPHERE_MEAN_DRIFT_mass_properties.csv",
                          delimiter=",")
        return hydro, body, wave, coefficients, mass
    with np.load(PORTABLE / "source.npz") as source:
        arrays = tuple(source[key].copy() for key in
                       ("body", "wave", "coefficients", "mass"))
    return (PORTABLE / "sphere.h5", *arrays)


def _ramp(time):
    return np.where(time < 20, (1 - np.cos(np.pi * time / 20)) / 2, 1)


def _radiation_force(hydro_file, velocity, output_axis):
    """Reconstruct the finite-window radiation integral from source HDF5."""
    dt = 0.01
    with h5py.File(hydro_file) as source:
        path = source["body1/hydro_coeffs/radiation_damping/impulse_response_fun"]
        lag = np.asarray(path["t"]).ravel()
        kernel = 1000 * CubicSpline(lag, np.asarray(path["K"]), axis=2)(
            np.arange(1001) * dt,
        )
    result = np.zeros(len(velocity))
    for input_axis in (0, 2, 4):
        tap = kernel[output_axis, input_axis]
        speed = velocity[:, input_axis]
        result += dt * fftconvolve(speed, tap, mode="full")[:len(speed)]
        result -= dt / 2 * tap[0] * speed
        result[1000:] -= dt / 2 * tap[-1] * speed[:-1000]
    return result


def test_sphere_control_surface_drift_matches_matlab_motion_and_forces():
    hydro, source, source_wave, coefficients, mass = _paired_data()
    assert (source.shape, source_wave.shape, coefficients.shape, mass.shape) == (
        (10_001, 43), (10_001, 2), (6, 3), (4,),
    )
    assert all(np.isfinite(values).all() for values in
               (source, source_wave, coefficients, mass))
    np.testing.assert_allclose(source_wave[:, 0], source[:, 0], rtol=0, atol=1e-8)
    wec = WEC("Sphere mean drift")
    sphere = wec.body("sphere", hydro, inertia=mass[1:],
                      mean_drift="control_surface")
    for axis in ("surge", "heave", "pitch"):
        wec.coordinate(axis, sphere.move(axis))
    result = wec.run(
        RegularCICWave(0.1, 2), dt=0.01, end_time=100,
        ramp_time=20, radiation_memory=10,
    )
    response = result.raw
    assert result.case["bodies"][0]["mean_drift"] == "control_surface"
    assert response.body_position.shape == (10_001, 1, 6)
    np.testing.assert_allclose(response.time, source[:, 0], rtol=0, atol=1e-8)
    np.testing.assert_allclose(response.wave_elevation, source_wave[:, 1],
                               rtol=0, atol=1e-12)

    for axis, label, position_limit, velocity_limit in (
        (0, "surge", 1e-4, 1e-4),
        (2, "heave", 7e-4, 2e-3),
        (4, "pitch", 1e-6, 1e-6),
    ):
        position_error = np.max(np.abs(
            response.body_position[:, 0, axis] - source[:, 1 + axis],
        ))
        velocity_error = np.max(np.abs(
            response.body_velocity[:, 0, axis] - source[:, 7 + axis],
        ))
        assert position_error < position_limit, f"{label} position {position_error}"
        assert velocity_error < velocity_limit, f"{label} velocity {velocity_error}"
    np.testing.assert_allclose(response.body_position[:, 0, [1, 3, 5]],
                               source[:, [2, 4, 6]], rtol=0, atol=1e-10)
    np.testing.assert_allclose(response.body_velocity[:, 0, [1, 3, 5]],
                               source[:, [8, 10, 12]], rtol=0, atol=1e-10)

    amplitude = 0.05
    ramp = _ramp(response.time)
    phase = np.pi * response.time
    expected_excitation = ramp[:, None] * (
        amplitude * (np.cos(phase[:, None]) * coefficients[:, 0]
                     - np.sin(phase[:, None]) * coefficients[:, 1])
        + amplitude**2 * coefficients[:, 2]
    )
    np.testing.assert_allclose(expected_excitation, source[:, 19:25],
                               rtol=0, atol=2e-9)
    outputs = dict(response.extra_outputs)
    np.testing.assert_allclose(outputs["body1_excitation_force"],
                               source[:, 19:25], rtol=0, atol=2e-9)
    drift = outputs["body1_mean_drift_force"]
    np.testing.assert_allclose(drift, amplitude**2 * ramp[:, None]
                               * coefficients[:, 2], rtol=0, atol=2e-9)
    for axis in (0, 2, 4):
        reconstructed = _radiation_force(
            hydro, response.body_velocity[:, 0], axis,
        )
        error = np.max(np.abs(reconstructed - source[:, 25 + axis]))
        assert error < {0: 0.25, 2: 5.0, 4: 0.001}[axis], (
            f"{('surge', 'heave', 'pitch')[(0, 2, 4).index(axis)]} "
            f"radiation force {error}"
        )


def test_sphere_mean_drift_cli_labels_each_force_component(tmp_path):
    hydro = _paired_data()[0]
    wec = WEC("Sphere mean drift")
    sphere = wec.body("sphere", hydro, inertia=(837.75804096,) * 3,
                      mean_drift="control_surface")
    for axis in ("surge", "heave", "pitch"):
        wec.coordinate(axis, sphere.move(axis))
    case = wec.to_case(RegularCICWave(0.1, 2), dt=0.01,
                       end_time=0, ramp_time=20, radiation_memory=10)
    case_file = tmp_path / "case.json"
    case_file.write_text(json.dumps(case), encoding="utf-8")
    output = tmp_path / "motion.csv"
    cli_main([str(case_file), "--output", str(output)])
    header = output.read_text(encoding="utf-8").splitlines()[0].split(",")
    values = np.loadtxt(output, delimiter=",", skiprows=1)
    assert len(header) == values.size
    for kind in ("mean_drift_force", "excitation_force"):
        for axis in ("surge", "sway", "heave", "roll", "pitch", "yaw"):
            assert f"body1_{kind}_{axis}" in header

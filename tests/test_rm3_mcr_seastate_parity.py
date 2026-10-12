"""Pair the three published imported-spectrum RM3 MCR cases with MATLAB."""

import os
from pathlib import Path

import h5py
import numpy as np
import pytest
from scipy.interpolate import CubicSpline
from scipy.signal import fftconvolve

from wecsim import (ImportedSpectrumWave, WEC, mcr_spectrum_files,
                    run_rm3_spectrum_mcr)
from wecsim.environment.irregularWave import imported_spectrum_components, synthesize_irregular_response


APPLICATIONS = os.environ.get("WEC_SIM_APPLICATIONS_DIR")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
pytestmark = pytest.mark.skipif(
    not (APPLICATIONS and REFERENCE),
    reason="paired MATLAB RM3 imported-spectrum output not provided",
)


def _max_error(actual, expected, limit, label):
    assert actual.shape == expected.shape, label
    assert np.isfinite(expected).all(), label
    error = float(np.max(np.abs(actual - expected)))
    assert error < limit, f"{label}: max error {error:.6g} exceeds {limit}"


def test_matlab_applied_added_mass_matches_pinned_coefficients():
    """Pair the matrix used by Simscape with the BEM matrix and mass split."""
    root = Path(APPLICATIONS)
    reference = Path(REFERENCE)
    hydro = root / "_Common_Input_Files/RM3/hydroData/rm3.h5"
    input_inertias = ((20_907_301, 21_306_090.66, 37_085_481.11),
                     (94_419_614.57, 94_407_091.24, 28_542_224.82))
    with h5py.File(hydro) as h5:
        rho = float(np.asarray(h5["simulation_parameters/rho"]).item())
        for body in (1, 2):
            prefix = f"body{body}"
            original = rho * np.asarray(
                h5[f"{prefix}/hydro_coeffs/added_mass/inf_freq"]
            )[:, 6 * (body - 1):6 * body]
            applied = original.copy()
            mass_shift = 2 * np.trace(original[:3, :3])
            applied[:3, :3] -= np.eye(3) * mass_shift
            applied[3:6, 3:6] -= np.diag(np.diag(original[3:6, 3:6]))
            for row, column in ((3, 4), (3, 5), (4, 5)):
                applied[row, column] -= original[row, column]
                applied[column, row] -= original[row, column]
            saved_original = np.loadtxt(reference / (
                f"RM3_MCR_SEASTATE_body{body}_added_mass_original.csv"),
                delimiter=",")
            saved_applied = np.loadtxt(reference / (
                f"RM3_MCR_SEASTATE_body{body}_added_mass_applied.csv"),
                delimiter=",")
            np.testing.assert_allclose(saved_original, original,
                                       rtol=1e-12, atol=1e-6)
            np.testing.assert_allclose(saved_applied, applied,
                                       rtol=1e-12, atol=1e-6)
            volume = float(np.asarray(h5[f"{prefix}/properties/disp_vol"]).item())
            nominal_mass = rho * volume
            properties = np.loadtxt(reference / (
                f"RM3_MCR_SEASTATE_body{body}_mass_properties.csv"),
                delimiter=",")
            expected_properties = np.r_[
                nominal_mass, nominal_mass, nominal_mass + mass_shift,
                input_inertias[body - 1],
                np.asarray(input_inertias[body - 1]) + np.diag(original)[3:6],
            ]
            np.testing.assert_allclose(properties, expected_properties,
                                       rtol=1e-12, atol=1e-6)
            nominal_rigid = np.diag(np.r_[
                np.repeat(properties[0], 3), properties[3:6],
            ])
            adjusted_rigid = np.diag(np.r_[
                np.repeat(properties[2], 3), properties[6:9],
            ])
            active = np.ix_([0, 2, 4], [0, 2, 4])
            np.testing.assert_allclose(
                (adjusted_rigid + saved_applied)[active],
                (nominal_rigid + saved_original)[active],
                rtol=1e-12, atol=1e-6,
            )


def test_matlab_added_mass_uses_extrapolated_acceleration():
    """Identify the acceleration reaching Simulink's tiny Transport Delay.

    Its output at each major time step extrapolates the two preceding
    accelerations to t - 1e-7 s. MATLAB postprocessing adds the pitch inertia
    removed from the applied matrix back to its reported added-mass force.
    The HDF5 matrix is not symmetric, so acceleration multiplies its transpose
    when the samples are represented as rows here.
    """
    reference = Path(REFERENCE)
    active = np.array([0, 2, 4])
    dt = 0.1
    delay = 1e-7  # The pinned rigid-body block sets 10e-8 s.
    extrapolation = 1 - delay / dt
    for body in (1, 2):
        original = np.loadtxt(
            reference / f"RM3_MCR_SEASTATE_body{body}_added_mass_original.csv",
            delimiter=",",
        )[np.ix_(active, active)]
        applied = np.loadtxt(
            reference / f"RM3_MCR_SEASTATE_body{body}_added_mass_applied.csv",
            delimiter=",",
        )[np.ix_(active, active)]
        for case in (1, 2, 3):
            trace = np.loadtxt(
                reference / f"RM3_MCR_SEASTATE_case{case}_body{body}_forces.csv",
                delimiter=",",
            )
            acceleration = trace[:, 19:25][:, active]
            delayed = (acceleration[1:-1] + extrapolation
                       * (acceleration[1:-1] - acceleration[:-2]))
            expected = delayed @ applied.T
            expected[:, 2] += acceleration[2:, 2] * original[2, 2]
            _max_error(
                expected, trace[2:, 7:13][:, active], 1e-6,
                f"case {case} body {body} delayed added-mass force",
            )


def test_matlab_sea_state_joint_and_force_balance():
    """Identify the joint motion and applied mass convention in the baseline.

    MATLAB moves twice the translational added-mass trace into each Simscape
    body. Its exported forceTotal retains the matching modified force, so the
    adjusted mass is needed to close the logged translational balance.
    """
    root = Path(APPLICATIONS)
    reference = Path(REFERENCE)
    hydro = root / "_Common_Input_Files/RM3/hydroData/rm3.h5"
    adjusted_mass = []
    centers = []
    inertias = (21_306_090.66, 94_407_091.24)
    with h5py.File(hydro) as h5:
        rho = float(np.asarray(h5["simulation_parameters/rho"]).item())
        for body in (1, 2):
            prefix = f"body{body}"
            volume = float(np.asarray(h5[f"{prefix}/properties/disp_vol"]).item())
            centers.append(float(np.asarray(h5[f"{prefix}/properties/cg"]).ravel()[2]))
            added = rho * np.asarray(
                h5[f"{prefix}/hydro_coeffs/added_mass/inf_freq"]
            )[:, 6 * (body - 1):6 * body]
            adjusted_mass.append(rho * volume + 2 * np.trace(added[:3, :3]))

    for case in (1, 2, 3):
        bodies = [np.loadtxt(reference / f"RM3_MCR_SEASTATE_case{case}_body{body}.csv",
                             delimiter=",") for body in (1, 2)]
        forces = [np.loadtxt(reference / f"RM3_MCR_SEASTATE_case{case}_body{body}_forces.csv",
                             delimiter=",") for body in (1, 2)]
        pto = np.loadtxt(reference / f"RM3_MCR_SEASTATE_case{case}_pto1.csv",
                         delimiter=",")
        angle = bodies[0][:, 5]
        sine, cosine = np.sin(angle), np.cos(angle)
        separation = centers[0] - centers[1] + pto[:, 3]
        _max_error(bodies[0][:, 1] - bodies[1][:, 1],
                   separation * sine, 1e-9, f"case {case} joint x")
        _max_error(bodies[0][:, 3] - bodies[1][:, 3],
                   separation * cosine, 1e-9, f"case {case} joint z")
        _max_error(pto[:, 15], -1_200_000 * pto[:, 9],
                   1e-6, f"case {case} PTO damper")
        _max_error(pto[:, 21], pto[:, 15] * pto[:, 9],
                   1e-6, f"case {case} PTO mechanical power")

        residual = np.zeros((len(angle), 4))
        for index in (0, 1):
            radius = bodies[index][:, 3] / cosine
            imbalance_x = (adjusted_mass[index] * forces[index][:, 19]
                           - bodies[index][:, 13])
            imbalance_z = (adjusted_mass[index] * forces[index][:, 21]
                           - bodies[index][:, 15])
            imbalance_pitch = (inertias[index] * forces[index][:, 23]
                               - bodies[index][:, 17])
            residual[:, 0] += imbalance_x
            residual[:, index + 1] += sine * imbalance_x + cosine * imbalance_z
            residual[:, 3] += (radius * cosine * imbalance_x
                               - radius * sine * imbalance_z
                               + imbalance_pitch)
        residual[:, 1] -= pto[:, 15]
        residual[:, 2] += pto[:, 15]
        _max_error(residual, np.zeros_like(residual), 1e-3,
                   f"case {case} generalized force balance")


def test_all_three_imported_waves_and_excitation_against_matlab():
    root = Path(APPLICATIONS)
    reference = Path(REFERENCE)
    mat_file = root / "Multiple_Condition_Runs/RM3_MCROPT3_SeaState/mcrExample.mat"
    files = mcr_spectrum_files(mat_file)
    assert tuple(path.name for path in files) == tuple(
        f"spectrumData{index}.mat" for index in (1, 2, 3)
    )
    hydro = root / "_Common_Input_Files/RM3/hydroData/rm3.h5"
    for index, spectrum_file in enumerate(files, start=1):
        components = imported_spectrum_components(hydro, spectrum_file)
        expected_components = np.loadtxt(
            reference / f"RM3_MCR_SEASTATE_case{index}_components.csv",
            delimiter=",",
        )
        assert expected_components.shape == (len(components.omega), 4)
        actual_components = np.column_stack((
            components.omega, components.spectral_amplitude,
            components.d_omega, components.phase[:, 0],
        ))
        np.testing.assert_allclose(actual_components, expected_components,
                                   rtol=0, atol=1e-11)
        assert len(components.directions) == 1
        assert components.directions[0] == 0
        expected_wave = np.loadtxt(
            reference / f"RM3_MCR_SEASTATE_case{index}_wave.csv",
            delimiter=",",
        )
        assert expected_wave.shape == (4001, 2)
        for body in (1, 2):
            expected = np.loadtxt(
                reference / f"RM3_MCR_SEASTATE_case{index}_body{body}.csv",
                delimiter=",",
            )
            assert expected.shape == (4001, 25)
            incident = synthesize_irregular_response(
                hydro, components, dt=0.1, end_time=400,
                ramp_time=100, body_number=body,
            )
            np.testing.assert_allclose(incident.time, expected_wave[:, 0],
                                       rtol=0, atol=1e-10)
            if body == 1:
                _max_error(incident.elevation, expected_wave[:, 1],
                           1e-10, f"case {index} wave elevation")
            for dof, name, force_limit in (
                (0, "surge", 1), (2, "heave", 1), (4, "pitch", 1),
            ):
                _max_error(incident.excitation_force[:, dof],
                           expected[:, 19 + dof], force_limit,
                           f"case {index} body {body} {name} excitation")


def test_public_python_wec_replays_published_imported_seas():
    """Exercise the Python builder against the three saved MATLAB sea states.

    This checks imported-wave propagation through the general WEC runner;
    the one-body heave motion here is not the published RM3 four-DOF joint.
    """
    root = Path(APPLICATIONS)
    reference = Path(REFERENCE)
    hydro = "_Common_Input_Files/RM3/hydroData/rm3.h5"
    inputs = root / "Multiple_Condition_Runs/RM3_MCROPT3_SeaState"
    for index in (1, 2, 3):
        wec = WEC(f"Imported RM3 sea {index}")
        body = wec.body("float", hydro)
        wec.coordinate("heave", body.move("heave"))
        wave = ImportedSpectrumWave(
            f"Multiple_Condition_Runs/RM3_MCROPT3_SeaState/"
            f"spectrumData{index}.mat"
        )
        result = wec.run(
            wave, dt=0.1, end_time=400, ramp_time=100,
            radiation_memory=15, base_dir=root,
        )
        source = np.loadtxt(
            reference / f"RM3_MCR_SEASTATE_case{index}_wave.csv",
            delimiter=",",
        )[:len(result.time)]
        _max_error(result.time, source[:, 0], 1e-10,
                   f"case {index} public WEC time")
        _max_error(result.wave_elevation, source[:, 1], 1e-10,
                   f"case {index} public WEC imported elevation")
        assert result.raw.auxiliary_files == (
            (inputs / f"spectrumData{index}.mat").resolve(),
        )


def test_matlab_radiation_force_on_matlab_velocities():
    """Pair the source convolution force independently of Python dynamics."""
    root = Path(APPLICATIONS)
    reference = Path(REFERENCE)
    hydro = root / "_Common_Input_Files/RM3/hydroData/rm3.h5"
    with h5py.File(hydro) as h5:
        rho = float(np.asarray(h5["simulation_parameters/rho"]).item())
        for body in (1, 2):
            irf = h5[f"body{body}/hydro_coeffs/radiation_damping/impulse_response_fun"]
            irf_time = np.asarray(irf["t"]).ravel()
            own = slice(6 * (body - 1), 6 * body)
            raw_kernel = np.asarray(irf["K"])[:, own, :]
            lag_count = 601  # Published 60 s memory at 0.1 s output steps.
            dt = 0.1
            kernel = rho * CubicSpline(irf_time, raw_kernel, axis=2)(
                np.arange(lag_count) * dt,
            )
            for case in (1, 2, 3):
                name = f"RM3_MCR_SEASTATE_case{case}_body{body}"
                motion = np.loadtxt(reference / f"{name}.csv", delimiter=",")
                force = np.loadtxt(reference / f"{name}_forces.csv", delimiter=",")
                velocity = motion[:, 7:13]
                calculated = np.zeros_like(velocity)
                for output in range(6):
                    for input_dof in range(6):
                        calculated[:, output] += dt * fftconvolve(
                            velocity[:, input_dof], kernel[output, input_dof],
                            mode="full",
                        )[:len(velocity)]
                calculated -= dt / 2 * np.einsum(
                    "ij,tj->ti", kernel[:, :, 0], velocity,
                )
                calculated[lag_count - 1:] -= dt / 2 * np.einsum(
                    "ij,tj->ti", kernel[:, :, -1],
                    velocity[:-(lag_count - 1)],
                )
                _max_error(calculated, force[:, 1:7], 1e-6,
                           f"case {case} body {body} radiation force")


def test_matlab_restoring_and_logged_force_decomposition():
    """Exclude unreported hydrodynamic terms from the sea-state mismatch."""
    root = Path(APPLICATIONS)
    reference = Path(REFERENCE)
    hydro = root / "_Common_Input_Files/RM3/hydroData/rm3.h5"
    with h5py.File(hydro) as h5:
        rho = float(np.asarray(h5["simulation_parameters/rho"]).item())
        g = float(np.asarray(h5["simulation_parameters/g"]).item())
        for body in (1, 2):
            prefix = f"body{body}"
            restoring = rho * g * np.asarray(
                h5[f"{prefix}/hydro_coeffs/linear_restoring_stiffness"]
            )
            center = np.asarray(h5[f"{prefix}/properties/cg"]).ravel()
            for case in (1, 2, 3):
                name = f"RM3_MCR_SEASTATE_case{case}_body{body}"
                motion = np.loadtxt(reference / f"{name}.csv", delimiter=",")
                force = np.loadtxt(reference / f"{name}_forces.csv", delimiter=",")
                displacement = motion[:, 1:7].copy()
                displacement[:, :3] -= center
                active = [0, 2, 4]  # RM3 floating-joint surge, heave, pitch.
                _max_error((displacement @ restoring.T)[:, active],
                           force[:, 13:19][:, active], 1e-6,
                           f"case {case} body {body} active restoring")
                total = (motion[:, 19:25] - force[:, 1:7]
                         - force[:, 7:13] - force[:, 13:19])
                _max_error(total, motion[:, 13:19], 1e-6,
                           f"case {case} body {body} total force")


def test_published_three_sea_state_mcr_against_matlab():
    root = Path(APPLICATIONS)
    reference = Path(REFERENCE)
    mat_file = root / "Multiple_Condition_Runs/RM3_MCROPT3_SeaState/mcrExample.mat"
    files = mcr_spectrum_files(mat_file)
    hydro = root / "_Common_Input_Files/RM3/hydroData/rm3.h5"
    result = run_rm3_spectrum_mcr(
        hydro, mat_file, added_mass_scheme="simulink_delay",
    )
    assert result.spectrum_files == files
    assert len(result.traces) == 3
    summary = np.loadtxt(reference / "RM3_MCR_SEASTATE_summary.csv", delimiter=",")
    assert summary.shape == (3,)

    for index, _ in enumerate(files, start=1):
        trace = result.traces[index - 1]
        response = trace.response
        assert response is not None
        assert response.body_position.shape == (4001, 2, 6)
        for body in (1, 2):
            expected = np.loadtxt(
                reference / f"RM3_MCR_SEASTATE_case{index}_body{body}.csv",
                delimiter=",",
            )
            assert expected.shape == (4001, 25)
            np.testing.assert_allclose(response.time, expected[:, 0],
                                       rtol=0, atol=1e-10)
            for dof, name, position_limit, velocity_limit in (
                (0, "surge", 0.020, 0.005),
                (2, "heave", 0.0015, 0.0013),
                (4, "pitch", 0.0009, 0.0004),
            ):
                _max_error(response.body_position[:, body - 1, dof],
                           expected[:, 1 + dof], position_limit,
                           f"case {index} body {body} {name} position")
                _max_error(response.body_velocity[:, body - 1, dof],
                           expected[:, 7 + dof], velocity_limit,
                           f"case {index} body {body} {name} velocity")

        pto = np.loadtxt(
            reference / f"RM3_MCR_SEASTATE_case{index}_pto1.csv",
            delimiter=",",
        )
        assert pto.shape == (4001, 25)
        np.testing.assert_allclose(response.time, pto[:, 0],
                                   rtol=0, atol=1e-10)
        center_gap = (response.body_position[0, 0, 2]
                      - response.body_position[0, 1, 2])
        pitch = response.body_position[:, 0, 4]
        cosine = np.cos(pitch)
        stroke = ((response.body_position[:, 0, 2]
                   - response.body_position[:, 1, 2]) / cosine - center_gap)
        speed = ((response.body_velocity[:, 0, 2]
                  - response.body_velocity[:, 1, 2]
                  + (center_gap + stroke) * np.sin(pitch)
                  * response.body_velocity[:, 0, 4]) / cosine)
        np.testing.assert_allclose(stroke, response.pto_stroke, rtol=0, atol=1e-12)
        np.testing.assert_allclose(speed, response.pto_velocity, rtol=0, atol=1e-12)
        _max_error(stroke, pto[:, 3], 0.0035, f"case {index} PTO stroke")
        _max_error(speed, pto[:, 9], 0.002,
                   f"case {index} PTO speed")
        _max_error(response.pto_force, pto[:, 15], 2_500,
                   f"case {index} PTO force")
        _max_error(trace.absorbed_power, -pto[:, 21], 1_800,
                   f"case {index} PTO absorbed power")
        np.testing.assert_allclose(-summary[index - 1],
                                   np.mean(-pto[1999:, 21]), rtol=0, atol=1e-6)

    _max_error(result.mean_absorbed_power, -summary, 75,
               "three sea-state MCR mean powers")

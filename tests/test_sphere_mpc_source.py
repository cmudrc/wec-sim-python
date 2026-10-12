"""Pair the published Sphere MPC waves and prediction matrices with MATLAB."""

import os
from pathlib import Path

import h5py
import numpy as np
import pytest
from scipy.io import loadmat

from wecsim import JONSWAPWave, run_sphere_mpc
from wecsim.environment.irregularWave import (
    jonswap_equal_energy_components, synthesize_irregular_response,
)
from wecsim.control.mpcControl import (
    integrate_sphere_mpc_force, predict_sphere_excitation,
    solve_sphere_mpc_qp,
)
from wecsim.control.mpcPlant import (
    build_sphere_mpc_matrices, replay_sphere_mpc_plant,
)
from wecsim.environment.waveClass import WaveClass


APPLICATIONS = os.environ.get("WEC_SIM_APPLICATIONS_DIR")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
pytestmark = pytest.mark.skipif(
    not (APPLICATIONS and REFERENCE),
    reason="published MATLAB Sphere MPC output and HDF5 absent",
)


def _load(name):
    return np.loadtxt(Path(REFERENCE) / f"SPHERE_MPC_{name}.csv",
                      delimiter=",", ndmin=2)


def _match(actual, expected, limit, label):
    assert actual.shape == expected.shape, label
    assert np.isfinite(actual).all() and np.isfinite(expected).all(), label
    error = float(np.max(np.abs(actual - expected)))
    assert error < limit, f"{label}: {error:.6g} exceeds {limit}"


def test_published_sphere_mpc_waves_and_prediction_model():
    apps = Path(APPLICATIONS)
    hydro = apps / "_Common_Input_Files/Sphere/hydroData/sphere.h5"
    coefficients = apps / "Controls/MPC/coeff.mat"
    source = loadmat(Path(REFERENCE) / "SPHERE_MPC_setup.mat",
                     simplify_cells=True)
    matrices = build_sphere_mpc_matrices(hydro, coefficients)
    for name in ("A", "Bu", "Bv", "C"):
        _match(getattr(matrices, name), np.asarray(source["plant"][name]),
               1e-10, f"prediction plant {name}")
    for name in ("Sx", "Su", "Sv", "Q", "H"):
        _match(getattr(matrices, name), np.asarray(source["setup"][name]),
               1e-10, f"MPC {name}")
    assert matrices.Sx.shape == (90, 9)
    assert matrices.Su.shape == (90, 31)
    assert np.linalg.eigvalsh((matrices.H + matrices.H.T) / 2).min() > 0

    components = _load("components")
    assert components.shape == (500, 4)
    phase = components[:, 3, None]
    assert JONSWAPWave(2.5, 8, phase_file="phases.csv").as_case()["type"] == "jonswap"
    sea = jonswap_equal_energy_components(
        hydro, significant_height=2.5, peak_period=8,
        directions=np.array([0]), spreading=np.array([1]), phase=phase,
    )
    _match(sea.omega, components[:, 0], 1e-12, "public JONSWAP frequencies")
    _match(sea.spectral_amplitude, components[:, 1], 1e-12,
           "public JONSWAP spectrum")
    _match(sea.d_omega, components[:, 2], 1e-12, "public JONSWAP widths")
    with h5py.File(hydro) as h5:
        frequencies = np.asarray(h5["/simulation_parameters/w"]).ravel()
    wave = WaveClass("irregular")
    wave.T = 8
    wave.H = 2.5
    wave.spectrumType = "JS"
    wave.freqDisc = "EqualEnergy"
    wave.numFreq = 500
    wave.phaseData = phase
    wave.waveDir = [0]
    wave.waveSpread = [1]
    wave.waveSetup([frequencies.min(), frequencies.max()], "infinite",
                   100, .01, 40000, 9.81, 1000, 400)
    assert wave.gamma == 1.0
    _match(np.asarray(wave.w).ravel(), components[:, 0], 1e-12,
           "JONSWAP frequencies")
    _match(np.asarray(wave.A).ravel(), components[:, 1], 1e-12,
           "JONSWAP amplitudes")
    _match(np.asarray(wave.dw).ravel(), components[:, 2], 1e-12,
           "JONSWAP widths")
    _match(np.asarray(wave.waveAmpTime).T, _load("wave"), 3e-12,
           "JONSWAP elevation")

    controller = _load("controller")
    body = _load("MPC_body1")
    incident = synthesize_irregular_response(
        hydro, sea, dt=.01, end_time=400, ramp_time=100,
    )
    _match(incident.elevation, _load("wave")[:, 1], 3e-12,
           "public JONSWAP elevation")
    _match(incident.excitation_force[:, 2], body[:, 21], 1e-7,
           "public Sphere excitation force")
    plant = _load("plant_output")
    assert controller.shape == (40001, 13)
    assert body.shape == (40001, 25)
    assert plant.shape == (40001, 4)
    _match(plant[:, 0], body[:, 0], 1e-12, "plant sample time")
    _match(plant[:, 3], controller[:, 3], 1e-8,
           "internal plant PTO force")
    replay = replay_sphere_mpc_plant(
        matrices, controller[:, 3], body[:, 21], dt=.01,
    )
    _match(replay[:, 0], plant[:, 1], 7e-5, "internal plant heave speed")
    _match(replay[:, 1], plant[:, 2], 5e-5, "internal plant heave position")
    _match(controller[:, 9], controller[:, 3] * body[:, 9], 1e-7,
           "logged controller power")
    assert controller[np.flatnonzero(np.abs(controller[:, 3]) > 1e-9)[0], 0] == 205.51
    assert 2e6 < np.max(np.abs(controller[:, 3])) < 2.5e6
    assert np.max(np.abs(np.diff(controller[:, 3]) / .01)) < 1.500001e6

    state = _load("full_state")
    forecast = _load("excitation_prediction")
    rate = _load("command_rate")
    iteration = _load("iteration")
    assert state.shape == (40001, 10)
    assert forecast.shape == (801, 32)
    assert rate.shape == iteration.shape == (801, 2)
    _match(state[:, [1, 2, 9]], plant[:, 1:], 1e-8,
           "logged full prediction-plant state")
    assert rate[np.flatnonzero(np.abs(rate[:, 1]) > 1e-9)[0], 0] == 205.0
    assert iteration[0, 1] == iteration[1, 1] == 1
    assert iteration[-1, 1] == 401

    for index in range(204, len(forecast)):
        sample = index * 50  # 0.5 s MPC sample on a 0.01 s body grid
        history = body[sample - 204 * 50:sample + 1:50, 21]
        prediction = predict_sphere_excitation(history)
        _match(prediction, forecast[index, 1:], 1e-6,
               f"excitation forecast at {forecast[index, 0]:.1f} s")

    for index in range(410, len(rate)):
        sample = index * 50
        plan, feasible = solve_sphere_mpc_qp(
            matrices, state[sample, 1:], forecast[index, 1:],
        )
        assert feasible, f"source MPC QP infeasible at {rate[index, 0]:.1f} s"
        _match(plan[:1], rate[index, 1:2], 2.0,
               f"force-rate command at {rate[index, 0]:.1f} s")

    reconstructed_force = integrate_sphere_mpc_force(
        rate[:, 1], dt=.01, end_time=400,
    )
    _match(reconstructed_force, controller[:, 3], 1e-6,
           "delayed and integrated MPC force")

    independent = run_sphere_mpc(
        hydro, coefficients, phase=phase,
    )
    assert independent.feasible.all()
    _match(independent.time, body[:, 0], 1e-12, "closed-loop time")
    _match(independent.wave_elevation, _load("wave")[:, 1], 3e-12,
           "closed-loop elevation")
    _match(independent.excitation_force, body[:, 21], 1e-7,
           "closed-loop heave excitation")
    _match(independent.command_rate, rate[:, 1], 25,
           "closed-loop MPC command")
    _match(independent.pto_force, controller[:, 3], 25,
           "closed-loop PTO force")
    _match(independent.position, body[:, 3], .002,
           "closed-loop physical heave")
    _match(independent.velocity, body[:, 9], .002,
           "closed-loop physical speed")
    _match(-independent.absorbed_power, controller[:, 9], 5000,
           "closed-loop source-signed power")
    _match(independent.internal_state[:, 0], plant[:, 1], 7e-5,
           "closed-loop internal speed")
    _match(independent.internal_state[:, 1], plant[:, 2], 5e-5,
           "closed-loop internal position")
    # MATLAB applies its first PTO force at 205.51 s, after the controller
    # has commanded a force rate at 205.00 s and passed a 0.5 s transition.
    assert np.count_nonzero(independent.pto_force[:20551]) == 0

    limited = run_sphere_mpc(
        hydro, coefficients, phase=phase, end_time=210,
        max_force_rate=400_000,
    )
    assert limited.feasible.all()
    assert np.max(np.abs(limited.command_rate)) <= 400_000.001
    assert independent.command_rate[410] < -400_000
    assert abs(limited.position[-1] - independent.position[21000]) > .1

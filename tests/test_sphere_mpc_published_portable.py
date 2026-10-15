"""Gate the published 400 s Sphere MPC closed loop in ordinary Python CI."""

from hashlib import sha256
from pathlib import Path

import numpy as np
from scipy.io import loadmat

from wecsim import run_sphere_mpc
from wecsim.control.mpcControl import predict_sphere_excitation
from wecsim.control.mpcPlant import build_sphere_mpc_matrices
from wecsim.environment.irregularWave import (
    jonswap_equal_energy_components, synthesize_irregular_response,
)


FIXTURE = Path(__file__).parent / "fixtures/matlab/sphere_mpc"
HYDRO = FIXTURE / "sphere.h5"
COEFFICIENTS = FIXTURE / "coeff.mat"
SETUP = FIXTURE / "setup.mat"
SOURCE = FIXTURE / "source.npz"
HASHES = {
    HYDRO: "0893381b0083306cbcbcbb8a4a6aef43818de107b7dd2229f638828510bee804",
    COEFFICIENTS: "98e470dc9f480456f5d774ae2ec3ba9593eddc4bd991dc0d13f93ffb8c7a2350",
    SETUP: "60d87933ee7245d013786ad1fae05228976f821fdbee4d9b6f8b005f4a5f12ab",
    SOURCE: "77a04a4b35b159ba9332216ce5f0b1c9a551d84401a95d447bbdd6fe34e8d0a3",
}


def _within(actual, expected, limit, label):
    actual = np.asarray(actual)
    expected = np.asarray(expected)
    assert actual.shape == expected.shape, label
    assert np.isfinite(actual).all() and np.isfinite(expected).all(), label
    error = float(np.max(np.abs(actual - expected)))
    assert error < limit, f"{label}: {error:.6g} exceeds {limit:.6g}"


def test_published_sphere_mpc_wave_forecast_and_closed_loop():
    for path, digest in HASHES.items():
        assert sha256(path.read_bytes()).hexdigest() == digest
    with np.load(SOURCE, allow_pickle=False) as archive:
        source = {name: archive[name] for name in archive.files}
    np.testing.assert_allclose(source["time"], np.arange(40001) * .01,
                               rtol=0, atol=1e-10)
    assert source["components"].shape == (500, 4)
    assert source["forecast"].shape == (801, 31)
    assert source["command_rate"].shape == (801,)

    matrices = build_sphere_mpc_matrices(HYDRO, COEFFICIENTS)
    setup = loadmat(SETUP, simplify_cells=True)
    for name in ("A", "Bu", "Bv", "C"):
        _within(getattr(matrices, name), np.asarray(setup["plant"][name]),
                1e-10, f"published plant matrix {name}")
    for name in ("Sx", "Su", "Sv", "Q", "H"):
        _within(getattr(matrices, name), np.asarray(setup["setup"][name]),
                1e-10, f"published MPC matrix {name}")

    components = source["components"]
    phase = components[:, 3, None]
    sea = jonswap_equal_energy_components(
        HYDRO, significant_height=2.5, peak_period=8,
        directions=(0,), spreading=(1,), phase=phase,
    )
    for actual, expected in (
        (sea.omega, components[:, 0]),
        (sea.spectral_amplitude, components[:, 1]),
        (sea.d_omega, components[:, 2]),
    ):
        _within(actual, expected, 1e-12, "published JONSWAP components")
    incident = synthesize_irregular_response(
        HYDRO, sea, dt=.01, end_time=400, ramp_time=100,
    )
    _within(incident.elevation, source["wave"], 3e-12,
            "published wave elevation")
    _within(incident.excitation_force[:, 2], source["excitation"],
            1e-7, "published heave excitation")

    # The source MPC samples the excitation every 0.5 s. Rebuild its 40-step
    # forecast from that saved incident force before closing the Python loop.
    for index in range(204, len(source["forecast"])):
        sample = index * 50
        history = source["excitation"][sample - 204 * 50:sample + 1:50]
        prediction = predict_sphere_excitation(history)
        _within(prediction, source["forecast"][index], 1e-6,
                f"excitation forecast at {index * .5:g} s")

    result = run_sphere_mpc(HYDRO, COEFFICIENTS, phase=phase)
    assert result.feasible.all()
    _within(result.time, source["time"], 1e-12, "time")
    _within(result.wave_elevation, source["wave"], 3e-12, "closed-loop wave")
    _within(result.excitation_force, source["excitation"],
            1e-7, "closed-loop excitation")
    _within(result.command_rate, source["command_rate"],
            25, "MPC command rate")
    _within(result.pto_force, source["pto_force"], 25, "PTO force")
    _within(result.position, source["position"], .002, "sphere heave")
    _within(result.velocity, source["velocity"], .002, "sphere heave speed")
    _within(-result.absorbed_power, source["pto_power"],
            5000, "source-signed PTO power")
    _within(result.internal_state[:, 0], source["internal_speed"],
            7e-5, "MPC plant speed")
    _within(result.internal_state[:, 1], source["internal_position"],
            5e-5, "MPC plant position")

    source_work = np.trapezoid(source["pto_power"], x=source["time"])
    python_work = np.trapezoid(-result.absorbed_power, x=result.time)
    assert source_work < 0
    assert abs(python_work - source_work) / abs(source_work) < .001

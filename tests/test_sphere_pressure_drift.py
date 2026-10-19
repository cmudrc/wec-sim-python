"""Pair WEC-Sim's pressure-integration spectral mean-drift option."""

from pathlib import Path

import h5py
import numpy as np
from scipy.interpolate import CubicSpline

from wecsim import PMWave, WEC
from wecsim.environment.irregularWave import pm_equal_energy_components


FIXTURE = Path(__file__).parent / "fixtures/matlab/sphere_pressure_drift"
HYDRO = Path(__file__).parent / "fixtures/matlab/sphere_mean_drift/sphere.h5"


def test_pressure_integration_drift_against_pinned_matlab():
    with np.load(FIXTURE / "source.npz") as source:
        body, wave, recorded, mass = (source[key] for key in
                                      ("body", "wave", "components", "mass"))
    assert (body.shape, wave.shape, recorded.shape, mass.shape) == (
        (3001, 19), (3001, 2), (500, 4), (4,),
    )
    assert all(np.isfinite(values).all() for values in
               (body, wave, recorded, mass))

    components = pm_equal_energy_components(
        HYDRO, significant_height=.5, peak_period=8,
        directions=[0], spreading=[1], seed=1, phase_generator="matlab",
    )
    for actual, expected in (
        (components.omega, recorded[:, 0]),
        (components.spectral_amplitude, recorded[:, 1]),
        (components.d_omega, recorded[:, 2]),
        (components.phase[:, 0], recorded[:, 3]),
    ):
        np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-12)

    # Compute the selected mean-drift term directly from the HDF5 rather
    # than asking the production irregular-wave builder for its result.
    with h5py.File(HYDRO) as h5:
        frequency = np.asarray(h5["simulation_parameters/w"]).ravel()
        coefficient = np.asarray(h5[
            "body1/hydro_coeffs/mean_drift/pressure_integration/val"
        ])[:, 0, :]
    expected_drift = (
        CubicSpline(frequency, coefficient, axis=1)(components.omega)
        @ (components.spectral_amplitude * components.d_omega)
    ) * 1000 * 9.81
    assert np.max(abs(expected_drift)) > 1

    wec = WEC("Sphere pressure-integration mean drift")
    sphere = wec.body("sphere", HYDRO, inertia=mass[1:],
                      mean_drift="pressure_integration")
    for axis in ("surge", "heave", "pitch"):
        wec.coordinate(axis, sphere.move(axis))
    result = wec.run(
        PMWave(.5, 8, seed=1, phase_generator="matlab"),
        dt=.01, end_time=30, ramp_time=10, radiation_memory=10,
    )
    output = dict(result.raw.extra_outputs)
    ramp = np.where(result.time < 10,
                    (1 - np.cos(np.pi * result.time / 10)) / 2, 1)
    np.testing.assert_allclose(result.time, body[:, 0], rtol=0, atol=1e-9)
    np.testing.assert_allclose(result.wave_elevation, wave[:, 1],
                               rtol=0, atol=1e-10)
    np.testing.assert_allclose(output["body1_mean_drift_force"],
                               ramp[:, None] * expected_drift,
                               rtol=0, atol=1e-6)
    assert np.max(abs(output["body1_excitation_force"]
                      - body[:, 13:19])) < 1e-5
    for axis, position_limit, speed_limit in (
        (0, .003, .004),
        (2, .003, .004),
        (4, 1e-4, 5e-5),
    ):
        position = result.bodies["sphere"].position[:, axis]
        speed = result.bodies["sphere"].velocity[:, axis]
        assert np.max(abs(position - body[:, 1 + axis])) < position_limit
        assert np.max(abs(speed - body[:, 7 + axis])) < speed_limit

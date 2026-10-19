"""Pair spectral mean drift on a derived Sphere Mean_Drift PM sea."""

from pathlib import Path

import h5py
import numpy as np
import pytest
from scipy.interpolate import CubicSpline

from wecsim import PMWave, WEC
from wecsim.environment.irregularWave import pm_equal_energy_components


FIXTURES = Path(__file__).parent / "fixtures/matlab"
HYDRO = Path(__file__).parent / "fixtures/matlab/sphere_mean_drift/sphere.h5"


def _device(mean_drift):
    wec = WEC("Sphere irregular spectral mean drift")
    sphere = wec.body("sphere", HYDRO, inertia=(837.75804096,) * 3,
                      mean_drift=mean_drift)
    for axis in ("surge", "heave", "pitch"):
        wec.coordinate(axis, sphere.move(axis))
    return wec


@pytest.mark.parametrize(
    ("option", "fixture"),
    (("control_surface", "sphere_irregular_mean_drift"),
     ("momentum_conservation", "sphere_irregular_momentum_drift")),
)
def test_irregular_spectral_mean_drift_against_pinned_matlab(option, fixture):
    with np.load(FIXTURES / fixture / "source.npz") as source:
        body = source["body"]
        wave = source["wave"]
        recorded_components = source["components"]
        mass = source["mass"]
    assert body.shape == (3001, 19)
    assert wave.shape == (3001, 2)
    assert recorded_components.shape == (500, 4)
    assert mass.shape == (4,)
    assert all(np.isfinite(values).all() for values in
               (body, wave, recorded_components, mass))
    np.testing.assert_allclose(mass[1:], 837.75804096, rtol=0, atol=1e-6)
    if option == "momentum_conservation":
        # The two independent MATLAB exports differ only in drift option;
        # pin their complete sea and rigid-body inputs before motion pairing.
        with np.load(FIXTURES / "sphere_irregular_mean_drift/source.npz") as control:
            np.testing.assert_array_equal(wave, control["wave"])
            np.testing.assert_array_equal(recorded_components,
                                          control["components"])
            np.testing.assert_array_equal(mass, control["mass"])

    components = pm_equal_energy_components(
        HYDRO, significant_height=0.5, peak_period=8,
        directions=[0], spreading=[1], seed=1,
        phase_generator="matlab",
    )
    for actual, expected in (
        (components.omega, recorded_components[:, 0]),
        (components.spectral_amplitude, recorded_components[:, 1]),
        (components.d_omega, recorded_components[:, 2]),
        (components.phase[:, 0], recorded_components[:, 3]),
    ):
        np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-12)

    # This is the pinned irregExcF mean-drift term, evaluated independently
    # from the force builder on the selected HDF5 drift coefficients.
    with h5py.File(HYDRO) as h5:
        frequency = np.asarray(h5["simulation_parameters/w"]).ravel()
        drift = np.asarray(h5[
            f"body1/hydro_coeffs/mean_drift/{option}/val"
        ])[:, 0, :]
    coefficient = CubicSpline(frequency, drift, axis=1)(components.omega)
    expected_drift = (coefficient @ (
        components.spectral_amplitude * components.d_omega
    )) * 1000 * 9.81
    # Momentum conservation has a smaller (but still material) surge term;
    # retain the existing control-surface floor for its original fixture.
    assert np.max(np.abs(expected_drift)) > (
        1 if option == "control_surface" else .5
    )

    result = _device(option).run(
        PMWave(0.5, 8, seed=1, phase_generator="matlab"),
        dt=.01, end_time=30, ramp_time=10, radiation_memory=10,
    )
    outputs = dict(result.raw.extra_outputs)
    ramp = np.where(result.time < 10,
                    (1 - np.cos(np.pi * result.time / 10)) / 2, 1)
    np.testing.assert_allclose(result.time, body[:, 0], rtol=0, atol=1e-9)
    np.testing.assert_allclose(result.wave_elevation, wave[:, 1],
                               rtol=0, atol=1e-10)
    np.testing.assert_allclose(outputs["body1_mean_drift_force"],
                               ramp[:, None] * expected_drift,
                               rtol=0, atol=1e-6)
    assert np.max(np.abs(outputs["body1_excitation_force"]
                         - body[:, 13:19])) < 1e-5

    for axis, position_limit, speed_limit in (
        (0, .003, .004),
        (2, .003, .004),
        (4, 1e-4, 5e-5),
    ):
        position = result.bodies["sphere"].position[:, axis]
        speed = result.bodies["sphere"].velocity[:, axis]
        assert np.max(np.abs(position - body[:, 1 + axis])) < position_limit
        assert np.max(np.abs(speed - body[:, 7 + axis])) < speed_limit

    without_drift = _device("none").run(
        PMWave(0.5, 8, seed=1, phase_generator="matlab"),
        dt=.01, end_time=30, ramp_time=10, radiation_memory=10,
    )
    assert abs(result.bodies["sphere"].position[-1, 0]
               - without_drift.bodies["sphere"].position[-1, 0]) > .05
    if option == "control_surface":
        assert abs(result.bodies["sphere"].position[-1, 4]
                   - without_drift.bodies["sphere"].position[-1, 4]) > .02
    else:
        # This HDF5 momentum coefficient acts in surge only; a pitch-effect
        # floor from the control-surface fixture would be inappropriate.
        np.testing.assert_allclose(expected_drift[1:], 0, rtol=0, atol=1e-12)

"""Small, deterministic checks of the broadband yaw force preparation."""

import numpy as np
import pytest

from wecsim.environment.irregularWave import IrregularComponents
from wecsim.hydrodynamics.passiveYaw import SampledPassiveYawExcitation


def _inputs():
    headings = np.array([0., 90., 180., 270.])
    frequency = np.array([1., 3.])
    real = np.empty((6, 4, 2))
    imaginary = np.empty_like(real)
    for dof in range(6):
        for heading in range(4):
            real[dof, heading] = (dof + 1) * (heading + 1) * np.array([1., 2.])
            imaginary[dof, heading] = (dof + 1) * (heading + 1) * np.array([.2, .6])
    hydro = {
        "simulation_parameters": {"wave_dir": headings, "w": frequency},
        "hydro_coeffs": {"excitation": {"re": real, "im": imaginary}},
    }
    components = IrregularComponents(
        omega=np.array([1.5, 2.5]),
        spectral_amplitude=np.array([.3, .5]),
        d_omega=np.array([.2, .4]),
        directions=np.array([0., 45.]),
        spreading=np.array([.25, .75]),
        phase=np.array([[.1, .4], [.7, 1.1]]),
    )
    return hydro, components


def test_broadband_yaw_force_sums_seas_ramps_and_wraps_headings():
    hydro, components = _inputs()
    model = SampledPassiveYawExcitation.from_hydro_data(
        hydro, components, dt=.25, end_time=.5, ramp_time=.5,
        rho=1000, g=9.81,
    )

    # Calculate one grid sample directly from the two incident seas and two
    # frequencies, independently of the model's batched matrix products.
    at_time = .25
    ramp = .5
    expected_elevation = 0.
    expected_grid = np.zeros((2, 4, 6))
    for sea in range(2):
        for frequency_index, omega in enumerate(components.omega):
            amplitude = np.sqrt(components.spectral_amplitude[frequency_index]
                                * components.d_omega[frequency_index]
                                * components.spreading[sea])
            angle = omega * at_time + components.phase[frequency_index, sea]
            expected_elevation += ramp * amplitude * np.cos(angle)
            for heading in range(4):
                for dof in range(6):
                    real = np.interp(omega, hydro["simulation_parameters"]["w"],
                                     hydro["hydro_coeffs"]["excitation"]["re"][dof, heading])
                    imaginary = np.interp(
                        omega, hydro["simulation_parameters"]["w"],
                        hydro["hydro_coeffs"]["excitation"]["im"][dof, heading],
                    )
                    expected_grid[sea, heading, dof] += (
                        ramp * amplitude * 1000 * 9.81
                        * (real * np.cos(angle) - imaginary * np.sin(angle))
                    )
    np.testing.assert_allclose(model.elevation[1], expected_elevation,
                               rtol=0, atol=1e-14)
    np.testing.assert_allclose(model.force_grid[1], expected_grid,
                               rtol=1e-14, atol=1e-10)
    np.testing.assert_array_equal(model.force_grid[0], 0)

    # At 45 degrees yaw, the first sea wraps to 315 degrees and the second
    # samples the 0-degree BEM heading. Rotate both force pairs to world axes.
    local = .5 * (expected_grid[0, 3] + expected_grid[0, 0]) + expected_grid[1, 0]
    expected_world = local.copy()
    c = np.sqrt(.5)
    for first in (0, 3):
        expected_world[first] = c * (local[first] - local[first + 1])
        expected_world[first + 1] = c * (local[first] + local[first + 1])
    np.testing.assert_allclose(model.force(at_time, np.pi / 4), expected_world,
                               rtol=1e-14, atol=1e-10)
    with pytest.raises(ValueError, match="grid time"):
        model.force(.3, 0.)


def test_broadband_yaw_rejects_nonphysical_directional_weights():
    hydro, components = _inputs()
    invalid = IrregularComponents(
        components.omega, components.spectral_amplitude, components.d_omega,
        components.directions, np.array([1.1, -.1]), components.phase,
    )
    with pytest.raises(ValueError, match="components are inconsistent"):
        SampledPassiveYawExcitation.from_hydro_data(
            hydro, invalid, dt=.25, end_time=.5, ramp_time=.5,
            rho=1000, g=9.81,
        )

"""Analytic check of the RM3 fitted-radiation source convention."""

import h5py
import numpy as np
import pytest

from wecsim.hydrodynamics.radiation import (
    RadiationStateSpace, replay_rm3_fitted_radiation,
)


def _first_order_fit(path):
    """Two self channels and one cross-body channel with known step responses."""
    with h5py.File(path, "w") as h5:
        h5.create_dataset("simulation_parameters/rho", data=1025.0)
        for body in (1, 2):
            fit = h5.create_group(
                f"body{body}/hydro_coeffs/radiation_damping/state_space"
            )
            orders = np.zeros((6, 12), dtype=int)
            a = np.zeros((6, 12, 1, 1))
            b = np.zeros_like(a)
            c = np.zeros_like(a)
            channels = (
                ((0, 0, -2.0, 3.0, 5.0),
                 (0, 6, -4.0, 2.0, 7.0))
                if body == 1 else ((0, 6, -3.0, 4.0, 6.0),)
            )
            for output, source, pole, gain, scale in channels:
                orders[output, source] = 1
                a[output, source, 0, 0] = pole
                b[output, source, 0, 0] = gain
                c[output, source, 0, 0] = scale
            fit.create_dataset("it", data=orders)
            fit.create_dataset("A/all", data=a)
            fit.create_dataset("B/all", data=b)
            fit.create_dataset("C/all", data=c)
            # MATLAB's pinned body block omits D even when HDF5 stores it.
            fit.create_dataset("D/all", data=np.full((6, 12), 1e6))


@pytest.mark.parametrize("body_to_body", [False, True])
def test_fitted_radiation_replay_matches_first_order_step_response(
    tmp_path, body_to_body,
):
    hydro = tmp_path / "synthetic_rm3.h5"
    _first_order_fit(hydro)
    time = np.arange(21) * 0.05
    velocity = np.zeros((len(time), 2, 6))
    velocity[:, 0, 0] = 1.0
    velocity[:, 1, 0] = 2.0

    force = replay_rm3_fitted_radiation(
        hydro, time, velocity, body_to_body=body_to_body,
    )
    own_float = 1025 * 5 * 3 / 2 * (1 - np.exp(-2 * time))
    cross = 1025 * 7 * 2 * 2 / 4 * (1 - np.exp(-4 * time))
    own_spar = 1025 * 6 * 4 * 2 / 3 * (1 - np.exp(-3 * time))
    expected = np.zeros_like(force)
    expected[:, 0, 0] = own_float + (cross if body_to_body else 0)
    expected[:, 1, 0] = own_spar
    np.testing.assert_allclose(force, expected, rtol=1e-11, atol=1e-9)
    np.testing.assert_array_equal(force[0], 0)  # The HDF5 D term is ignored.

    # The reduced model must retain the same active-channel steady response.
    model = RadiationStateSpace.from_rm3_hdf5(
        hydro, body_to_body=body_to_body,
    )
    assert model.A.shape == ((3, 3) if body_to_body else (2, 2))
    speed = np.zeros(12)
    speed[[0, 6]] = [1, 2]
    steady = model.C @ np.linalg.solve(-model.A, model.B @ speed)
    np.testing.assert_allclose(
        steady[[0, 6]],
        [1025 * 5 * 3 / 2 + (1025 * 7 * 2 * 2 / 4 if body_to_body else 0),
         1025 * 6 * 4 * 2 / 3],
        rtol=1e-12, atol=1e-9,
    )

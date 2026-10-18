"""TurbSim BTS decoding and pinned MOST MATLAB reader comparison."""

import os
import struct

import numpy as np
import pytest
from scipy.io import loadmat

from wecsim import MostWindField, read_turbsim_bts
from wecsim.environment.turbSim import TurbSimWind


def test_small_bts_grid_order_scaling_and_tower(tmp_path):
    path = tmp_path / "small.bts"
    header = struct.pack(
        "<h4i12fi", 7, 2, 2, 1, 2, 5, 4, .5, 8, 20, 10,
        2, 10, 4, -2, 5, 1, 4,
    )
    # Within each time step: z, y, U/V/W, then one tower U/V/W.
    raw = np.arange(2 * (3 * 2 * 2 + 3), dtype="<i2")
    path.write_bytes(header + b"test" + raw.tobytes())
    wind = read_turbsim_bts(path)
    assert wind.velocity.shape == (2, 3, 2, 2)
    np.testing.assert_array_equal(wind.y, [-2, 2])
    np.testing.assert_array_equal(wind.z, [10, 15])
    np.testing.assert_allclose(wind.velocity[0, :, 0, 0], [-5, .75, .2])
    np.testing.assert_allclose(wind.velocity[0, :, 1, 1], [-.5, 3, 2])
    np.testing.assert_allclose(wind.tower_velocity[1, :, 0], [8.5, 7.5, 5.6])
    assert wind.dt == .5 and wind.hub_height == 20
    assert wind.description == "test"
    short = path.with_suffix(".short")
    short.write_bytes(path.read_bytes()[:-2])
    with pytest.raises(ValueError, match="file size"):
        read_turbsim_bts(short)


@pytest.mark.parametrize("dz, dy, dt", [
    (5., 4., 0.), (5., -4., .5), (np.nan, 4., .5),
])
def test_bts_rejects_invalid_spacing_before_advection(tmp_path, dz, dy, dt):
    path = tmp_path / "invalid-spacing.bts"
    header = struct.pack(
        "<h4i12fi", 7, 2, 2, 0, 2, dz, dy, dt, 8, 20, 10,
        2, 10, 4, -2, 5, 1, 0,
    )
    path.write_bytes(header + np.zeros(2 * 3 * 2 * 2, dtype="<i2").tobytes())
    with pytest.raises(ValueError, match="spatial or time metadata"):
        read_turbsim_bts(path)


def test_advection_rejects_zero_time_step_from_direct_wind_data():
    wind = TurbSimWind(np.zeros((30, 3, 2, 2)), np.empty((30, 3, 0)),
                       np.array([-100., 100.]), np.array([0., 1.]),
                       0., 10., 8., "")
    with pytest.raises(ValueError, match="time step"):
        MostWindField(wind)


def test_single_cell_bts_grid_allows_unused_zero_spacings(tmp_path):
    path = tmp_path / "single-cell.bts"
    header = struct.pack(
        "<h4i12fi", 7, 1, 1, 0, 2, 0, 0, .5, 8, 20, 10,
        2, 10, 4, -2, 5, 1, 0,
    )
    path.write_bytes(header + np.zeros(2 * 3, dtype="<i2").tobytes())
    wind = read_turbsim_bts(path)
    np.testing.assert_array_equal(wind.y, [0])
    np.testing.assert_array_equal(wind.z, [10])
    assert wind.dt == .5


@pytest.mark.skipif(
    not (os.environ.get("WEC_SIM_MOST_BTS") and os.environ.get("WEC_SIM_MOST_WIND_BASELINE")),
    reason="pinned MOST BTS and MATLAB reader output not provided",
)
def test_published_most_wind_against_matlab_reader():
    wind = read_turbsim_bts(os.environ["WEC_SIM_MOST_BTS"])
    reference = loadmat(os.environ["WEC_SIM_MOST_WIND_BASELINE"])
    nt, _, ny, nz = wind.velocity.shape
    assert (nt, ny, nz) == (20750, 12, 12)
    np.testing.assert_allclose(
        [nt, ny, nz, wind.dt, wind.hub_height],
        reference["metadata"].ravel(), rtol=0, atol=1e-12,
    )
    np.testing.assert_allclose(wind.y, reference["y"].ravel(), rtol=0, atol=1e-12)
    np.testing.assert_allclose(wind.z, reference["z"].ravel(), rtol=0, atol=1e-12)
    selected_y = [0, 5, 11]
    selected_z = [0, 6, 11]
    traces = wind.velocity[:, :, selected_y, :][:, :, :, selected_z]
    np.testing.assert_allclose(traces, reference["traces"], rtol=0, atol=1e-12)
    np.testing.assert_allclose(wind.velocity.sum(axis=(2, 3)),
                               reference["spatial_sum"], rtol=1e-13, atol=1e-10)
    np.testing.assert_allclose(np.square(wind.velocity).sum(axis=(2, 3)),
                               reference["spatial_square_sum"], rtol=1e-13, atol=1e-9)

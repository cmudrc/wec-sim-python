"""MOST's wind-plane advection against the pinned MATLAB preprocessor."""

import os
from pathlib import Path

import numpy as np
import pytest
from scipy.io import loadmat

from wecsim import MostWindField, read_turbsim_bts


def test_most_wind_advection_sample_order_and_bounds():
    # One component carries its time index so X-plane shifts are observable.
    from wecsim.environment.turbSim import TurbSimWind

    velocity = np.broadcast_to(np.arange(30)[:, None, None, None],
                               (30, 3, 2, 2)).copy()
    wind = TurbSimWind(velocity, np.empty((30, 3, 0)),
                       np.array([-100., 100.]), np.array([0., 1.]),
                       1., 10., 8., "")
    field = MostWindField(wind, speed=10)
    assert field.discarded == 10
    assert field.n_time == 11
    first = field.at_index(0)
    assert first.shape == (3, 10, 2, 2)
    np.testing.assert_array_equal(first[0, :, 0, 0], np.arange(13, 3, -1))
    sampler = field.sampler(0)
    np.testing.assert_allclose(
        sampler([-25, 0, .5]), first[:, :2].mean(axis=(1, 2, 3)),
        rtol=0, atol=1e-14,
    )
    np.testing.assert_allclose(
        field.sampler_at(.5)([-25, 0, .5]), sampler([-25, 0, .5]) + .5,
        rtol=0, atol=1e-14,
    )
    final = field.sampler(field.n_time - 1)([-25, 0, .5])
    for time in (field.time[-1], field.time[-1] + .1, 100.):
        np.testing.assert_array_equal(field.sampler_at(time)([-25, 0, .5]),
                                      final)
    for time in (-.1, np.nan, np.inf):
        with pytest.raises(ValueError, match="outside the available record"):
            field.sampler_at(time)
    np.testing.assert_array_equal(
        sampler([field.x[-1], wind.y[-1], wind.z[-1]]),
        first[:, -1, -1, -1],
    )
    assert np.isnan(sampler([field.x[-1]+1, 0, .5])).all()
    with pytest.raises(IndexError):
        field.at_index(field.n_time)


def test_single_cell_cross_section_interpolates_only_along_advection_axis():
    from wecsim.environment.turbSim import TurbSimWind

    components = np.array([1., 2., 3.])
    velocity = np.arange(30)[:, None, None, None] * components[None, :, None, None]
    wind = TurbSimWind(velocity, np.empty((30, 3, 0)),
                       np.array([100.]), np.array([2.]),
                       1., 10., 8., "")
    field = MostWindField(wind, speed=10)
    spatial_midpoint = field.at_index(0)[:, :2, 0, 0].mean(axis=1)
    np.testing.assert_allclose(field.sampler(0)([-25., 100., 2.]),
                               spatial_midpoint, rtol=0, atol=1e-14)
    np.testing.assert_allclose(field.sampler_at(.5)([-25., 100., 2.]),
                               spatial_midpoint + .5 * components,
                               rtol=0, atol=1e-14)


def test_wind_sampler_interpolates_all_three_spatial_axes_and_time():
    """An affine field has an exact trilinear and temporal interpolant.

    The three components use different gradients so swapping Y/Z, a vector
    component, or an advection plane changes the predicted blade wind.
    """
    from wecsim.environment.turbSim import TurbSimWind

    time = np.arange(30, dtype=float)[:, None, None, None]
    y = np.array([-100., 100.])
    z = np.array([0., 20.])
    time_gradient = np.array([1., -2., 3.])
    y_gradient = np.array([.1, .2, -.3])
    z_gradient = np.array([.4, -.5, .6])
    offset = np.array([8., 2., -1.])
    velocity = (
        time * time_gradient[None, :, None, None]
        + y[None, None, :, None] * y_gradient[None, :, None, None]
        + z[None, None, None, :] * z_gradient[None, :, None, None]
        + offset[None, :, None, None]
    )
    wind = TurbSimWind(velocity, np.empty((30, 3, 0)), y, z, 1., 10., 8., "")
    field = MostWindField(wind, speed=10)

    # At X=-30/-20 m, output t=2 selects TurbSim frames 15/14.
    # Halfway between those planes and a quarter-step later, the effective
    # source time is 14.75. The other coordinates are also interior points.
    point = [-25., 25., 7.]
    expected = (14.75 * time_gradient + 25. * y_gradient
                + 7. * z_gradient + offset)
    np.testing.assert_allclose(field.sampler_at(2.25)(point), expected,
                               rtol=0, atol=1e-13)
    np.testing.assert_allclose(field.sampler_at(2.)(point),
                               expected - .25 * time_gradient,
                               rtol=0, atol=1e-13)
    for outside in ([-25., 101., 7.], [-25., 25., 21.]):
        assert np.isnan(field.sampler_at(2.25)(outside)).all()


@pytest.mark.skipif(not os.environ.get("WEC_SIM_MOST_ADVECTION_DIR"),
                    reason="pinned MATLAB MOST advection output not provided")
def test_most_wind_planes_against_run_turbsim():
    directory = Path(os.environ["WEC_SIM_MOST_ADVECTION_DIR"])
    wind = read_turbsim_bts(directory / "WIND_8mps.bts")
    field = MostWindField(wind)
    reference = loadmat(directory.parent / "matlab-most-advection.mat")
    values = reference["values"]
    expected_times = int(os.environ.get("WEC_SIM_MOST_EXPECTED_WIND_TIMES", "249"))
    assert values.shape == (expected_times, 3, 10, 12, 12)
    assert field.n_time == expected_times and field.discarded == 376
    np.testing.assert_allclose(field.time, reference["time"].ravel(), rtol=0, atol=1e-12)
    np.testing.assert_array_equal(field.x, reference["x"].ravel())
    np.testing.assert_array_equal(wind.y, reference["y"].ravel())
    np.testing.assert_array_equal(wind.z, reference["z"].ravel())
    for index in range(field.n_time):
        np.testing.assert_allclose(field.at_index(index), values[index],
                                   rtol=0, atol=1e-12)

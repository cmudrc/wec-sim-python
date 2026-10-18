"""Physical wind vectors and domain boundaries used by MOST sampling."""

import numpy as np
import pytest

from wecsim import MostConstantWind, MostWindField


def test_constant_wind_direction_magnitude_and_domain():
    wind = MostConstantWind(12., end_time=5., direction=(3., 4., 0.),
                            domain_sizes=(20., 10., 8.))
    sample = wind.sampler_at(2.5)
    expected = [7.2, 9.6, 0.]
    np.testing.assert_allclose(sample([0., 0., 4.]), expected, rtol=0,
                               atol=1e-14)
    np.testing.assert_allclose(sample([10., -5., 8.]), expected, rtol=0,
                               atol=1e-14)
    # A blade element beyond any face of the wind box has no wind sample.
    for point in ([10.01, 0., 4.], [0., -5.01, 4.], [0., 0., -0.01],
                  [0., 0., 8.01]):
        assert np.isnan(sample(point)).all()
    with pytest.raises(ValueError, match="outside the available record"):
        wind.sampler_at(5.01)
    with pytest.raises(ValueError, match="finite 3-vector"):
        sample([np.inf, 0., 4.])


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

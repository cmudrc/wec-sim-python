"""Physical wind vectors and domain boundaries used by MOST sampling."""

import numpy as np
import pytest

from wecsim import MostConstantWind


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

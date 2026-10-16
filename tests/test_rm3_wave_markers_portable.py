"""Keep published RM3 marker elevations paired in ordinary Python CI."""

from hashlib import sha256
from pathlib import Path

import numpy as np

from wecsim import RegularWave


SOURCE = (Path(__file__).parent / "fixtures/matlab/rm3_wave_markers/source.npz")
SOURCE_SHA256 = "c767a06d2b641857dfe999cf32861a34509ba3a90fa5ced67e7f92848db90613"


def test_published_rm3_marker_elevations():
    """Rebuild all 25 world locations from the pinned wave input for 400 s."""
    assert sha256(SOURCE.read_bytes()).hexdigest() == SOURCE_SHA256
    with np.load(SOURCE, allow_pickle=False) as source:
        time = source["time"]
        locations = source["locations"]
        markers = source["markers"]
        origin = source["origin"]
        parameters = source["wave_parameters"]

    assert time.shape == (4001,)
    assert locations.shape == (25, 2)
    assert markers.shape == (4001, 25)
    assert origin.shape == (4001,)
    assert parameters.shape == (7,)
    assert all(np.isfinite(values).all() for values in
               (time, locations, markers, origin, parameters))
    np.testing.assert_allclose(time, np.arange(4001) * .1, rtol=0, atol=1e-10)
    np.testing.assert_array_equal(parameters[:6], [2.5, 8, 0, 100, 9.81, 200])

    grid = np.arange(-20, 21, 10)
    x, y = np.meshgrid(grid, grid)
    np.testing.assert_array_equal(
        locations, np.column_stack((x.ravel(order="F"), y.ravel(order="F"))),
    )
    height, period, direction, ramp, gravity, depth, source_k = parameters
    actual = RegularWave(height, period, direction).elevation_at(
        time, locations, water_depth=depth, ramp_time=ramp, g=gravity,
    )
    error = float(np.max(np.abs(actual - markers)))
    assert error < 1e-10, f"RM3 marker elevation: {error:.6g} m"
    center = np.flatnonzero(np.all(locations == [0, 0], axis=1))
    assert center.size == 1
    origin_error = float(np.max(np.abs(actual[:, center[0]] - origin)))
    assert origin_error < 1e-10, f"RM3 origin elevation: {origin_error:.6g} m"

    omega = 2 * np.pi / period
    dispersion_error = abs(
        omega**2 - gravity * source_k * np.tanh(source_k * depth)
    ) / omega**2
    assert dispersion_error < 1e-10

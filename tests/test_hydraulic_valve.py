"""Flow direction and continuity for the PTO-Sim four-port check valve."""

import numpy as np
import pytest

from wecsim.components.hydraulic import RectifyingCheckValve


def _valve():
    return RectifyingCheckValve(.61, .002, 1e-8, 1.5e6, 0,
                                850, 200, 8.137375096991404e-6)


def test_rectifying_valve_routes_high_cylinder_pressure():
    valve = _valve()
    a, b, high, low = valve.flows(20e6, 10e6, 15e6, 8e6)
    assert a < 0 and high > 0 and low < 0
    np.testing.assert_allclose(a + b + high + low, 0, atol=1e-16)
    assert valve.check_flow(5e6) > valve.check_flow(-5e6)
    assert valve.check_flow(0) == 0


def test_rectifying_valve_accepts_pressure_arrays():
    valve = _valve()
    flows = valve.flows(np.array([20e6, 10e6]), 10e6, 15e6, 8e6)
    assert all(flow.shape == (2,) for flow in flows)
    np.testing.assert_allclose(np.sum(flows, axis=0), 0, atol=1e-16)


def test_invalid_valve_geometry_is_rejected():
    with pytest.raises(ValueError, match="bounds"):
        RectifyingCheckValve(.61, .002, .003, 1.5e6, 0,
                             850, 200, 8e-6)

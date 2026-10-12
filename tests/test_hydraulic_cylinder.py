"""Pressure and force balances for the PTO-Sim compressible cylinder."""

import numpy as np
import pytest

from wecsim.components.hydraulic import CompressibleCylinder


def test_equal_chambers_exchange_pressure_under_piston_motion():
    cylinder = CompressibleCylinder(.0378, .0378, 1.86e9, 70, 35)
    rate_a, rate_b = cylinder.pressure_rates(0, .2, 0, 0)
    np.testing.assert_allclose(rate_a, 1.86e9 * .2 / 35)
    np.testing.assert_allclose(rate_b, -rate_a)
    np.testing.assert_allclose(cylinder.force(2e7, 1e7), -.0378 * 1e7)


def test_port_flow_can_balance_piston_displacement():
    cylinder = CompressibleCylinder(.0378, .0378, 1.86e9, 70, 35)
    rate_a, rate_b = cylinder.pressure_rates(
        np.array([0., 1.]), np.array([.1, -.2]),
        np.array([-.00378, .00756]),
        np.array([.00378, -.00756]),
    )
    np.testing.assert_allclose(rate_a, 0, atol=1e-8)
    np.testing.assert_allclose(rate_b, 0, atol=1e-8)


def test_invalid_cylinder_geometry_is_rejected():
    with pytest.raises(ValueError, match="positive"):
        CompressibleCylinder(0, .0378, 1.86e9, 70, 35)

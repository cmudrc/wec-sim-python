"""Geometry checks for the two published OSWEC motion converters."""

import numpy as np
import pytest

from wecsim.components.crank import AdjustableRodCrank, FixedRodCrank, PitchRodLinkage


@pytest.mark.parametrize("linkage", [
    FixedRodCrank(3, 1.3, 5),
    AdjustableRodCrank(3, 1.3, 5),
])
def test_stroke_starts_at_zero_and_jacobian_matches_finite_difference(linkage):
    assert abs(linkage.stroke_and_jacobian(0)[0]) < 1e-14
    angles = np.linspace(-.8, .8, 21)
    stroke, jacobian = linkage.stroke_and_jacobian(angles)
    delta = 1e-6
    plus = linkage.stroke_and_jacobian(angles + delta)[0]
    minus = linkage.stroke_and_jacobian(angles - delta)[0]
    np.testing.assert_allclose((plus - minus) / (2 * delta), jacobian,
                               rtol=0, atol=1e-9)
    assert np.isfinite(stroke).all()


def test_fixed_rod_rejects_unreachable_angle():
    linkage = FixedRodCrank(3, 1.3, 2)
    with pytest.raises(ValueError, match="outside"):
        linkage.stroke_and_jacobian(np.pi)


def test_body_local_pitch_rod_jacobian():
    linkage = PitchRodLinkage(
        anchor=(5.6021271782, -8.7), hinge=(0, -8.9),
        body_center=(0, -3.9), body_point=(.9, -3.1),
    )
    angles = np.linspace(-.7, .7, 31)
    stroke, jacobian = linkage.stroke_and_jacobian(angles)
    delta = 1e-6
    finite_difference = (
        linkage.stroke_and_jacobian(angles + delta)[0]
        - linkage.stroke_and_jacobian(angles - delta)[0]
    ) / (2 * delta)
    np.testing.assert_allclose(stroke[15], 0, rtol=0, atol=1e-14)
    np.testing.assert_allclose(jacobian, finite_difference,
                               rtol=0, atol=1e-9)

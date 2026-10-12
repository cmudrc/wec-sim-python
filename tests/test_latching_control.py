"""Check the pinned Sphere latching law before pairing body dynamics."""

import pytest

from wecsim import LatchingControl
from wecsim.control.controls import LatchingState


def test_reversal_switches_to_finite_dissipative_latch():
    control = LatchingControl(gain=49_181, latch_damping=37_308_296,
                              latch_time=2.4)
    state = LatchingState(elapsed_normal=0.21, previous_velocity=0.1)
    force, state = control.sample(-0.1, state, 0.01)
    assert force == pytest.approx(3_730_829.6)
    assert force * -0.1 < 0
    assert state.elapsed_latched == pytest.approx(0.01)
    assert state.elapsed_normal == 0
    for _ in range(239):
        force, state = control.sample(-0.1, state, 0.01)
        assert force == pytest.approx(3_730_829.6)
    assert state.elapsed_latched == pytest.approx(2.4)
    # Repeated floating-point addition can leave the timer a few ulps below
    # 2.4 s, just as in the pinned MATLAB function.
    force, state = control.sample(-0.1, state, 0.01)
    assert force == pytest.approx(3_730_829.6)
    force, state = control.sample(-0.1, state, 0.01)
    assert force == pytest.approx(4_918.1)
    assert state.elapsed_latched == 0
    assert state.elapsed_normal == pytest.approx(0.01)


def test_minimum_normal_interval_suppresses_retrigger():
    control = LatchingControl(100, 1_000, 0.5)
    state = LatchingState(elapsed_normal=0.2, previous_velocity=1)
    force, state = control.sample(-2, state, 0.01)
    assert force == 200
    assert state.elapsed_latched == 0


@pytest.mark.parametrize("normal,latch,duration", [
    (0, 100, 1), (100, 0, 1), (100, 100, 0),
    (float("nan"), 100, 1), (100, float("inf"), 1),
])
def test_reject_invalid_settings(normal, latch, duration):
    with pytest.raises(ValueError):
        LatchingControl(normal, latch, duration)

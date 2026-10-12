"""Verify the published discrete PTO switching law independently of motion."""

import pytest

from wecsim.control.controls import DeclutchingControl, DeclutchingState


def test_declutch_after_reversal_and_reengage_after_duration():
    control = DeclutchingControl(gain=232_020, declutch_time=0.8)
    state = DeclutchingState(elapsed_on=0.21, previous_velocity=0.01)
    force, state = control.sample(-0.01, state, 0.01)
    assert force == 0
    assert state.elapsed_off == pytest.approx(0.01)
    assert state.elapsed_on == 0
    for _ in range(79):
        force, state = control.sample(-0.02, state, 0.01)
        assert force == 0
    assert state.elapsed_off == pytest.approx(0.8)
    force, state = control.sample(-0.02, state, 0.01)
    assert force == pytest.approx(4_640.4)
    assert state.elapsed_off == 0
    assert state.elapsed_on == pytest.approx(0.01)


def test_normal_time_suppresses_immediate_retrigger():
    control = DeclutchingControl(gain=100, declutch_time=0.8)
    state = DeclutchingState(elapsed_on=0.2, previous_velocity=1)
    force, state = control.sample(-2, state, 0.01)
    assert force == 200
    assert state.elapsed_off == 0


@pytest.mark.parametrize("gain,duration", [(-1, 0.8), (0, 0.8),
                                          (100, 0), (float("nan"), 1)])
def test_invalid_declutching_settings(gain, duration):
    with pytest.raises(ValueError):
        DeclutchingControl(gain, duration)

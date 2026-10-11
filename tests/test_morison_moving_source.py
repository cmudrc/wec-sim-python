"""Pinned MATLAB regWaveMorison option-1 force at prescribed moving states."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim.morison import MorisonElement, regular_morison_source_force


REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")


@pytest.mark.skipif(not REFERENCE, reason="paired MATLAB Morison states absent")
def test_moving_regular_source_force():
    rows = np.loadtxt(Path(REFERENCE) / "MORISON_MOVING_REGULAR_states.csv",
                      delimiter=",")
    assert rows.shape == (8, 25)
    assert np.isfinite(rows).all()
    element = MorisonElement(
        point=(.8, -.5, 1.2),
        drag_coefficient=(1.1, .7, 1.3),
        added_mass_coefficient=(.8, 1.2, .6),
        area=(2, 3, 1.5), volume=1.7,
    )
    actual = np.stack([
        regular_morison_source_force(
            [element], time=row[0], position=row[1:7],
            velocity=row[7:13], acceleration=row[13:19],
            wave_height=2.5, wave_period=8, direction=25,
            water_depth=30, ramp_time=10, rho=1025,
        ) for row in rows
    ])
    error = np.max(np.abs(actual - rows[:, 19:25]))
    assert error < 1e-9, f"moving Morison source force differs by {error:.6g}"
    assert np.max(np.abs(actual[:-2])) > 1
    np.testing.assert_array_equal(actual[-2], 0)


@pytest.mark.skipif(not REFERENCE, reason="paired MATLAB Morison states absent")
def test_moving_regular_normal_tangential_source_force():
    rows = np.loadtxt(
        Path(REFERENCE) / "MORISON_MOVING_NORMAL_TANGENTIAL_states.csv",
        delimiter=",",
    )
    assert rows.shape == (8, 25)
    assert np.isfinite(rows).all()
    element = MorisonElement(
        point=(.8, -.5, 1.2),
        drag_coefficient=(1.1, .7, 1.3),
        added_mass_coefficient=(.8, 1.2, .6),
        area=(2, 3, 1.5), volume=1.7,
    )
    actual = np.stack([
        regular_morison_source_force(
            [element], time=row[0], position=row[1:7],
            velocity=row[7:13], acceleration=row[13:19],
            wave_height=2.5, wave_period=8, direction=25,
            water_depth=30, ramp_time=10, rho=1025,
            body_morison=2, element_axes=[(.3, -.7, 1.1)],
        ) for row in rows
    ])
    error = float(np.max(np.abs(actual - rows[:, 19:25])))
    assert error < 1e-8, f"normal/tangential source force differs by {error:.6g}"
    assert np.max(np.abs(actual[:-2])) > 1
    np.testing.assert_array_equal(actual[-2], 0)


@pytest.mark.skipif(not REFERENCE, reason="paired MATLAB Morison states absent")
@pytest.mark.parametrize("profile", ["uniform", "power", "linear"])
def test_moving_regular_source_current_force(profile):
    rows = np.loadtxt(
        Path(REFERENCE) / f"MORISON_MOVING_REGULAR_CURRENT_{profile}_states.csv",
        delimiter=",",
    )
    assert rows.shape == (8, 25)
    element = MorisonElement(
        point=(.8, -.5, 1.2), drag_coefficient=(1.1, .7, 1.3),
        added_mass_coefficient=(.8, 1.2, .6), area=(2, 3, 1.5), volume=1.7,
    )
    actual = np.stack([
        regular_morison_source_force(
            [element], time=row[0], position=row[1:7],
            velocity=row[7:13], acceleration=row[13:19],
            wave_height=2.5, wave_period=8, direction=25,
            water_depth=30, ramp_time=10, rho=1025,
            current_speed=.8, current_direction=45,
            current_profile=profile, current_depth=30,
        ) for row in rows
    ])
    error = float(np.max(np.abs(actual - rows[:, 19:25])))
    assert error < 1e-9, f"{profile} current force differs by {error:.6g}"
    assert np.max(np.abs(actual[:-2])) > 1
    np.testing.assert_array_equal(actual[-2], 0)

"""Pinned MATLAB regWaveMorison option-1 force at prescribed moving states."""

import os
from hashlib import sha256
from pathlib import Path

import numpy as np
import pytest

from wecsim.hydrodynamics.morison import (
    MorisonElement, regular_morison_source_force,
    regular_wave_normal_tangential_morison_terms,
)


REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
ARCHIVED = Path(__file__).parent / "fixtures/matlab/morison_moving_regular"
ARCHIVED_SHA256 = {
    "MORISON_MOVING_REGULAR_states.csv":
        "466ff5a912ca07062d699f9275e6ea2db2ae42b20ee3b6844a967e53f96ca827",
    "MORISON_MOVING_NORMAL_TANGENTIAL_states.csv":
        "ce0014df5dea6638b747b33c42166bfca6563c86e5fdc3c970de8f8221e0d81e",
    "MORISON_MOVING_REGULAR_CURRENT_uniform_states.csv":
        "a9bc61ea3822c9f39c7c97e5e08ebf4de3bb7c5abaec54dffc641f6f4328c48f",
    "MORISON_MOVING_REGULAR_CURRENT_power_states.csv":
        "3c69e1d2cdcf6effa396e5f0dc3938fc36c9749e5c769bb3b78ddc53ff0b30b2",
    "MORISON_MOVING_REGULAR_CURRENT_linear_states.csv":
        "8e9e01f85122e13c86c8176356d75aa8ac19feb7efaa382bace253b6fca3acff",
}


def _source_rows(filename):
    reference = Path(REFERENCE) / filename if REFERENCE else ARCHIVED / filename
    if not REFERENCE:
        assert sha256(reference.read_bytes()).hexdigest() == ARCHIVED_SHA256[filename]
    return np.loadtxt(reference, delimiter=",")


def test_moving_regular_source_force():
    rows = _source_rows("MORISON_MOVING_REGULAR_states.csv")
    assert rows.shape == (8, 25)
    assert np.isfinite(rows).all()
    assert np.any(np.linalg.norm(rows[:, 4:7], axis=1) > .1)
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


def test_moving_regular_normal_tangential_source_force():
    rows = _source_rows("MORISON_MOVING_NORMAL_TANGENTIAL_states.csv")
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


def test_physical_normal_tangential_force_on_unrotated_source_states():
    rows = _source_rows("MORISON_MOVING_NORMAL_TANGENTIAL_states.csv")
    element = MorisonElement(
        point=(.8, -.5, 1.2),
        drag_coefficient=(1.1, .7, 1.3),
        added_mass_coefficient=(.8, 1.2, .6),
        area=(2, 3, 1.5), volume=1.7,
        axis=(.3, -.7, 1.1),
    )
    for index in (0, 1, 6):
        row = rows[index]
        np.testing.assert_array_equal(row[4:7], 0)
        applied, mass = regular_wave_normal_tangential_morison_terms(
            [element], position=row[1:7], velocity=row[7:13],
            time=row[0], wave_height=2.5, wave_period=8,
            direction=25, water_depth=30, ramp_time=10, rho=1025,
        )
        np.testing.assert_allclose(
            applied - mass @ row[13:19], row[19:25],
            rtol=0, atol=1e-9,
        )
        np.testing.assert_allclose(mass, mass.T, rtol=0, atol=1e-10)
        assert np.linalg.eigvalsh(mass).min() > -1e-9


@pytest.mark.parametrize("profile", ["uniform", "power", "linear"])
def test_moving_regular_source_current_force(profile):
    rows = _source_rows(f"MORISON_MOVING_REGULAR_CURRENT_{profile}_states.csv")
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

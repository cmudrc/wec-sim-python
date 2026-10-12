"""Source option-1 current increment on a saved coupled Sphere trajectory."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim import Current, RegularWave
from wecsim.hydrodynamics.morison import MorisonElement, regular_morison_source_force


REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
pytestmark = pytest.mark.skipif(
    not REFERENCE, reason="paired moving-current MATLAB output absent",
)


def test_source_current_increment_on_coupled_motion():
    folder = Path(REFERENCE)
    body = np.loadtxt(
        folder / "SPHERE_MOVING_MORISON_CURRENT_1m-ME_body1.csv",
        delimiter=",",
    )
    delta = np.loadtxt(
        folder / "SPHERE_MOVING_MORISON_CURRENT_direct_delta.csv",
        delimiter=",",
    )
    wave = np.loadtxt(
        folder / "SPHERE_MOVING_MORISON_CURRENT_wave.csv",
        delimiter=",",
    )
    assert body.shape == (4001, 37)
    assert delta.shape == (4001, 7)
    assert wave.shape == (4001, 2)
    assert all(np.isfinite(values).all() for values in (body, delta, wave))
    np.testing.assert_allclose(delta[:, 0], body[:, 0], atol=1e-12, rtol=0)
    np.testing.assert_allclose(wave[:, 0], body[:, 0], atol=1e-12, rtol=0)
    predicted_wave = RegularWave(
        1, 8, current=Current(.8, 0, "power", 30),
    ).elevation_at(body[:, 0], [[0, 0]], ramp_time=10)[:, 0]
    assert np.max(np.abs(predicted_wave - wave[:, 1])) < 1e-10

    element = MorisonElement(
        point=(0, 0, -2), drag_coefficient=(0, 0, 1),
        added_mass_coefficient=(0, 0, 1), area=(0, 0, 100), volume=20,
    )
    reproduced = np.stack([
        regular_morison_source_force(
            [element], time=row[0], position=row[1:7],
            velocity=row[7:13], acceleration=row[31:37],
            wave_height=1, wave_period=8, direction=0,
            water_depth=1e6, ramp_time=10, rho=1000,
            current_speed=.8, current_direction=0,
            current_profile="power", current_depth=30,
        ) - regular_morison_source_force(
            [element], time=row[0], position=row[1:7],
            velocity=row[7:13], acceleration=row[31:37],
            wave_height=1, wave_period=8, direction=0,
            water_depth=1e6, ramp_time=10, rho=1000,
        )
        for row in body
    ])
    error = float(np.max(np.abs(reproduced - delta[:, 1:])))
    assert error < 1e-6, f"source current increment differs by {error:.6g}"
    assert np.max(np.abs(delta[:, 1])) > 1
    assert np.max(np.abs(delta[:, 5])) > 1

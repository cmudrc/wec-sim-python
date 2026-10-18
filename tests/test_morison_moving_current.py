"""Current loading on a tilted body-local axial element."""

from pathlib import Path

import numpy as np
import pytest

from wecsim import Current, RegularWave, WEC
from wecsim.hydrodynamics.morison import MorisonElement, regular_wave_axial_morison_terms


ELEMENT = MorisonElement(
    point=(0, 0, -2), drag_coefficient=(0, 0, 1),
    added_mass_coefficient=(0, 0, 0), area=(0, 0, 2), volume=1,
)


def _force(time, profile="uniform", depth=None, pitch=np.pi / 6):
    return regular_wave_axial_morison_terms(
        [ELEMENT], position=(0, 0, -2, 0, pitch, 0),
        velocity=np.zeros(6), time=time, wave_height=0, wave_period=8,
        ramp_time=10, water_depth=30, rho=1000,
        current_speed=1, current_direction=0,
        current_profile=profile, current_depth=depth,
    )


def test_tilted_axial_element_projects_horizontal_current():
    force, added = _force(10)
    axis = np.array([.5, 0, np.sqrt(3) / 2])
    np.testing.assert_allclose(force[:3], 250 * axis, atol=1e-12)
    np.testing.assert_allclose(force[3:], 0, atol=1e-12)
    np.testing.assert_array_equal(added, 0)
    np.testing.assert_allclose(_force(5)[0], force / 4, atol=1e-12)
    np.testing.assert_array_equal(_force(10, pitch=0)[0], 0)


def test_depth_profiles_and_public_wave_configuration():
    uniform = _force(10)[0]
    linear = _force(10, "linear", 30)[0]
    power = _force(10, "power", 30)[0]
    assert np.linalg.norm(uniform) > np.linalg.norm(power) > np.linalg.norm(linear)
    assert RegularWave(1, 8, current=Current(1, 30, "power", 30)).as_case()[
        "current"]["profile"] == "power"
    with pytest.raises(ValueError, match="current settings"):
        _force(10, "linear")


@pytest.mark.parametrize("current", [
    Current(True), Current(1, True), Current(1, 0, "power", True),
    Current(np.bool_(True)),
])
def test_public_wave_rejects_boolean_current_quantities(current):
    # A boolean must not become 1 m/s, 1 degree, or 1 m when the wave is run.
    with pytest.raises(ValueError, match="numeric speed, direction, and depth"):
        RegularWave(0, 8, current=current).as_case()


def test_public_runner_accepts_current_with_finite_depth_hydro():
    hydro = (Path(__file__).parent / "test_objects/test_bodyclass/testData"
             / "hydroData/ellipsoid.h5")
    wec = WEC("tilted ellipsoid")
    body = wec.body("ellipsoid", hydro, inertia=(1e6, 1e6, 1e6))
    for axis in ("surge", "heave", "pitch"):
        wec.coordinate(axis, body.move(axis))
    wec.morison_element(
        body, point=body.at(0, 0, -2), drag_coefficient=(0, 0, 1),
        added_mass_coefficient=(0, 0, 0), area=(0, 0, 2), volume=1,
    )
    run = dict(dt=.01, end_time=.02, ramp_time=0,
               initial_coordinate={"pitch": .1})
    still = wec.run(RegularWave(0, 8), **run)
    flowing = wec.run(RegularWave(0, 8, current=Current(.8, 0)), **run)
    assert np.max(np.abs(flowing.body_forces["ellipsoid"]
                         - still.body_forces["ellipsoid"])) > 1
    assert np.isfinite(flowing.bodies["ellipsoid"].position).all()

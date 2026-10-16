"""A configured OSWEC hinge must start at its supplied pitch state."""

from pathlib import Path

import numpy as np
import pytest

from wecsim import NoWave, PMWave, RegularWave, WEC, WorldPoint


ROOT = Path(__file__).resolve().parents[1]
HYDRO = ROOT / "examples/data/oswec.h5"
NONLINEAR_HYDRO = ROOT / "tests/fixtures/matlab/oswec_multi_wave/oswec.h5"
FLAP_MESH = ROOT / "tests/test_objects/test_bodyclass/testData/body_3_test/flap.stl"


def _hinge(*, nonlinear=False, source_solver=False):
    wec = WEC("OSWEC initial pitch")
    flap = wec.body(
        "flap", NONLINEAR_HYDRO if nonlinear else HYDRO,
        mass=127_000, inertia=(0, 1.85e6, 0),
        **({"geometry_file": FLAP_MESH, "nonlinear_hydro": "instantaneous"}
           if nonlinear else {}),
    )
    base = wec.fixed_body("base", center_gravity=(0, 0, -10.9), mass=999)
    wec.fixed_hinge(
        flap, base, location=WorldPoint(0, 0, -10),
        pto_location=WorldPoint(0, 0, -8.9), damping=12_000,
        stiffness=500,
        added_mass_scheme="simulink_delay" if source_solver else "implicit",
    )
    return wec


def test_regular_hinge_free_decay_matches_the_equivalent_named_coordinate():
    initial = dict(initial_coordinate={"pitch": .2},
                   initial_speed={"pitch": -.01})
    wave = RegularWave(0, 8)
    fixed = _hinge().run(wave, dt=.05, end_time=1, ramp_time=0, **initial)

    mapped = WEC("Same pitch with an explicit coordinate")
    flap = mapped.body("flap", HYDRO, mass=127_000, inertia=(0, 1.85e6, 0))
    pitch = mapped.coordinate(
        "pitch", flap.move("pitch", pivot=WorldPoint(0, 0, -8.9)),
    )
    mapped.rotational_pto("hinge", pitch, damping=12_000, stiffness=500)
    equivalent = mapped.run(wave, dt=.05, end_time=1, ramp_time=0, **initial)

    assert fixed.case["constraint"]["initial_coordinate"] == {"pitch": .2}
    assert fixed.case["constraint"]["initial_speed"] == {"pitch": -.01}
    np.testing.assert_array_equal(fixed.coordinates["pitch"].position,
                                  equivalent.coordinates["pitch"].position)
    np.testing.assert_array_equal(fixed.coordinates["pitch"].velocity,
                                  equivalent.coordinates["pitch"].velocity)
    np.testing.assert_array_equal(fixed.bodies["flap"].position,
                                  equivalent.bodies["flap"].position)
    np.testing.assert_array_equal(fixed.bodies["flap"].velocity,
                                  equivalent.bodies["flap"].velocity)
    np.testing.assert_array_equal(fixed.ptos["hinge"].force,
                                  equivalent.ptos["hinge"].force)
    assert fixed.coordinates["pitch"].position[0] == .2
    assert fixed.coordinates["pitch"].velocity[0] == -.01
    assert np.min(fixed.ptos["hinge"].absorbed_power) >= 0


def test_irregular_hinge_starts_memory_and_pto_from_supplied_state():
    wave = PMWave(1, 8, directions=(0,), spreading=(1,),
                  frequency_count=8, seed=7)
    initial = dict(initial_coordinate={"pitch": .12},
                   initial_speed={"pitch": -.02})
    result = _hinge().run(
        wave, dt=.05, end_time=.2, ramp_time=0,
        radiation_memory=.2, **initial,
    )
    zero = _hinge().run(
        wave, dt=.05, end_time=.2, ramp_time=0,
        radiation_memory=.2,
    )
    assert result.coordinates["pitch"].position[0] == .12
    assert result.coordinates["pitch"].velocity[0] == -.02
    np.testing.assert_allclose(result.ptos["hinge"].force[0],
                               -12_000 * -.02 - 500 * .12)
    assert np.max(np.abs(result.coordinates["pitch"].position
                         - zero.coordinates["pitch"].position)) > .1


def test_no_wave_hinge_free_decay_uses_radiation_and_initial_state():
    result = _hinge().run(
        NoWave(), dt=.05, end_time=.3, radiation_memory=.2,
        initial_coordinate={"pitch": .12}, initial_speed={"pitch": -.02},
    )
    angle = result.coordinates["pitch"].position
    speed = result.coordinates["pitch"].velocity
    assert result.case["wave"] == {"type": "none"}
    assert angle[0] == .12
    assert speed[0] == -.02
    assert np.isfinite(angle).all()
    assert np.max(np.abs(angle - angle[0])) > 0
    np.testing.assert_array_equal(result.wave_elevation, np.zeros(len(result.time)))
    np.testing.assert_allclose(result.ptos["hinge"].force,
                               -12_000 * speed - 500 * angle)

    # The independent coordinate-map route must produce the same free decay.
    mapped = WEC("Mapped still-water pitch")
    flap = mapped.body("flap", HYDRO, mass=127_000, inertia=(0, 1.85e6, 0))
    pitch = mapped.coordinate(
        "pitch", flap.move("pitch", pivot=WorldPoint(0, 0, -8.9)),
    )
    mapped.rotational_pto("hinge", pitch, damping=12_000, stiffness=500)
    equivalent = mapped.run(
        NoWave(), dt=.05, end_time=.3, radiation_memory=.2,
        initial_coordinate={"pitch": .12}, initial_speed={"pitch": -.02},
    )
    np.testing.assert_array_equal(angle, equivalent.coordinates["pitch"].position)
    np.testing.assert_array_equal(speed, equivalent.coordinates["pitch"].velocity)
    np.testing.assert_array_equal(result.ptos["hinge"].force,
                                  equivalent.ptos["hinge"].force)


@pytest.mark.parametrize("source_solver", (False, True))
def test_nonlinear_hinge_preserves_initial_pitch_and_speed(source_solver):
    result = _hinge(nonlinear=True, source_solver=source_solver).run(
        RegularWave(0, 8), dt=.1, end_time=.1, ramp_time=0,
        solver="wecsim_ode4" if source_solver else "implicit",
        initial_coordinate={"pitch": .08}, initial_speed={"pitch": -.015},
    )
    assert result.coordinates["pitch"].position[0] == .08
    assert result.coordinates["pitch"].velocity[0] == -.015
    np.testing.assert_allclose(result.ptos["hinge"].force[0],
                               -12_000 * -.015 - 500 * .08)
    assert np.isfinite(result.bodies["flap"].position).all()

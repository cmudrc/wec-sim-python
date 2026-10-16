"""Guard the physical parameter domains of paired PTO force elements."""

from dataclasses import replace

import numpy as np
import pytest

from wecsim import LinearHardStops, OrificePTO
from wecsim.dynamics.friction import StribeckFriction


@pytest.mark.parametrize("change", (
    {"lower_bound": .6},
    {"lower_stiffness": 0},
    {"upper_damping": -1},
    {"lower_transition_width": 0},
    {"upper_transition_width": float("nan")},
))
def test_hard_stops_reject_nonphysical_contact_parameters(change):
    stops = LinearHardStops(-.6, .6, 1e8, 1e8)
    with pytest.raises(ValueError, match="hard-stop"):
        replace(stops, **change)


def test_stribeck_friction_is_odd_and_dissipative_across_zero_speed():
    friction = StribeckFriction(
        breakaway_torque=12, breakaway_speed=.4,
        coulomb_torque=7, viscous_damping=2,
    )
    speeds = np.array([-100, -2, -.4, -.01, 0, .01, .4, 2, 100])
    torque = friction.torque(speeds)
    np.testing.assert_allclose(torque, -friction.torque(-speeds), atol=1e-14)
    assert np.all(torque * speeds <= 0)
    assert torque[4] == 0
    np.testing.assert_allclose(torque[-1], -207, atol=1e-12)


@pytest.mark.parametrize("change", (
    {"breakaway_torque": 6},
    {"breakaway_speed": 0},
    {"coulomb_torque": -1},
    {"viscous_damping": -1},
    {"breakaway_torque": float("inf")},
))
def test_stribeck_friction_rejects_nonpassive_parameters(change):
    friction = StribeckFriction(12, .4, 7, 2)
    with pytest.raises(ValueError, match="rotational friction"):
        replace(friction, **change)


def test_orifice_force_opposes_flow_and_flags_supersonic_use():
    orifice = OrificePTO(
        piston_area=.1, orifice_area=.01,
        discharge_coefficient=.5, air_density=1.2,
        mach_threshold=.3, sound_speed=100,
    )
    speeds = np.array([-4., -2., 0., 2., 4.])
    response = orifice.evaluate(speeds)
    np.testing.assert_array_equal(response.compressibility_flag,
                                  [True, False, False, False, True])
    assert np.all(response.force * speeds <= 0)
    np.testing.assert_allclose(response.pressure_drop[-1],
                               4 * response.pressure_drop[-2])
    np.testing.assert_allclose(response.absorbed_power,
                               -response.force * speeds)


@pytest.mark.parametrize("change", (
    {"orifice_area": 0},
    {"piston_area": True},
    {"air_density": float("nan")},
    {"mach_threshold": -1},
))
def test_orifice_rejects_nonphysical_parameters(change):
    orifice = OrificePTO(.1, .01)
    with pytest.raises(ValueError, match="orifice areas"):
        replace(orifice, **change)


def test_orifice_rejects_nonfinite_piston_speed():
    with pytest.raises(ValueError, match="piston speed must be finite"):
        OrificePTO(.1, .01).evaluate([0, float("nan")])

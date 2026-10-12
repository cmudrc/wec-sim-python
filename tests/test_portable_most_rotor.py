"""Small aerodynamic checks for the MOST rotor without external blade files."""

import numpy as np

from wecsim.components.mostRotor import MostRotor
from wecsim.control.mostController import MostBaselineController
from wecsim.environment.turbSim import MostConstantWind
from wecsim.hydrodynamics.mostBEM import MostBEM


def _small_rotor():
    # A smooth, symmetric airfoil keeps the BEM root solve physical while
    # letting the test distinguish rotor geometry from source table loading.
    angle = np.linspace(-180, 180, 721)
    foil = np.zeros((len(angle), 4, 1))
    foil[:, 0, 0] = angle
    foil[:, 1, 0] = .8 * np.sin(2 * np.deg2rad(angle))
    foil[:, 2, 0] = .05 + .1 * np.sin(np.deg2rad(angle))**2
    return MostBEM(
        radius=np.array([5., 10., 15.]), radius_interval=np.full(3, 5.),
        twist=np.zeros(3), chord=np.ones(3), curve_angle=np.zeros(3),
        sweep=np.zeros(3), curve=np.zeros(3),
        airfoil_index=np.zeros(3, dtype=int), airfoil=foil,
        hub_radius=3., tip_radius=17.,
    )


def test_three_blade_bem_rotates_loads_cyclically_and_gains_thrust_with_wind():
    bem = _small_rotor()
    hub = np.zeros(14)
    hub[13] = 1  # rotor speed, rad/s
    slow = bem.loads(hub, 0, [5, 0, 0])
    fast = bem.loads(hub, 0, [10, 0, 0])
    np.testing.assert_allclose(
        fast, bem.loads(hub, 0, lambda _point: np.array([10., 0., 0.])),
        rtol=0, atol=0,
    )
    assert np.isfinite(fast).all()
    assert np.sum(fast[0]) > np.sum(slow[0]) > 0
    np.testing.assert_array_equal(fast[2], 0)  # straight blade, no lateral force
    hub[12] = 2 * np.pi / 3
    np.testing.assert_allclose(bem.loads(hub, 0, [10, 0, 0]),
                               np.roll(fast, -1, axis=1), rtol=0, atol=1e-8)


def test_short_rotor_run_couples_bem_torque_and_baseline_control():
    rotor = MostRotor(
        bem=_small_rotor(), controller=MostBaselineController.iea15mw(),
        inertia=1e9, hub_offset=np.zeros(3),
        shaft_axis=np.array([1., 0., 0.]),
    )
    time = np.array([0., .05, .1])
    still = np.zeros((3, 6))
    result = rotor.simulate(time, still, still,
                            MostConstantWind(10, end_time=.1))
    assert result.blade_root_load.shape == (3, 6, 3)
    assert np.isfinite(result.blade_root_load).all()
    assert np.all(np.diff(result.azimuth) > 0)
    assert np.all(np.diff(result.rotor_speed) < 0)
    assert np.all(result.generator_torque >= 0)
    np.testing.assert_allclose(result.azimuth[1] - result.azimuth[0],
                               .05 * (result.rotor_speed[0]
                                      + result.rotor_speed[1]) / 2,
                               rtol=0, atol=1e-12)

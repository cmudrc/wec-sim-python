"""Small aerodynamic checks for the MOST rotor without external blade files."""

import numpy as np

from wecsim import MostCoupled, MostPlatformHydrodynamics, MostTowerReaction
from wecsim.components.mostRotor import MostRotor
from wecsim.control.mostController import MostBaselineController
from wecsim.dynamics.mostMooring import MostStaticMooring
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


def test_iea_blade_reader_preserves_geometry_and_distinct_airfoil_polars(tmp_path):
    """File conversion must preserve the local blade geometry used by BEM."""
    def write_table(name, rows, header_lines):
        np.savetxt(
            tmp_path / name, rows,
            header="\n".join(["source header"] * header_lines), comments="",
        )

    # Raw blade radii are relative to the 3.94 m hub. Airfoil indices in the
    # input are one-based, and adjacent stations use midpoint properties.
    write_table("IEA-15-240-RWT_blade.dat", np.array([
        [1, 0, 0, 0, 2, 1, 1],
        [3, .2, 1, 10, 4, 3, 2],
        [5, .4, 2, 20, 6, 5, 2],
    ]), 6)
    # The duplicate zero-angle row has deliberately different coefficients:
    # source preprocessing retains the first occurrence for interpolation.
    write_table("IEA-15-240-RWT_Airfoil_0.dat", np.array([
        [-20, -1, .1, -.2],
        [0, 0, .2, 0],
        [0, 99, 99, 99],
        [20, 1, .3, .2],
    ]), 54)
    write_table("IEA-15-240-RWT_Airfoil_1.dat", np.array([
        [-20, 0, 0, -.1],
        [0, .5, .05, 0],
        [20, 1, .1, .1],
    ]), 54)

    model = MostBEM.from_iea15mw(tmp_path)
    np.testing.assert_allclose(model.radius, [5.94, 7.94], rtol=0, atol=1e-12)
    np.testing.assert_allclose(model.radius_interval, [2, 2], rtol=0, atol=0)
    np.testing.assert_allclose([model.hub_radius, model.tip_radius],
                               [4.94, 8.94], rtol=0, atol=1e-12)
    np.testing.assert_allclose(model.twist, [3, 5], rtol=0, atol=0)
    np.testing.assert_allclose(model.chord, [2, 4], rtol=0, atol=0)
    np.testing.assert_allclose(model.curve_angle, np.deg2rad([5, 15]),
                               rtol=0, atol=1e-15)
    np.testing.assert_allclose(model.sweep, [.5, 1.5], rtol=0, atol=0)
    np.testing.assert_allclose(model.curve, [.1, .3], rtol=0, atol=1e-15)
    np.testing.assert_array_equal(model.airfoil_index, [0, 1])

    angles = np.linspace(-20, 20, 4)
    np.testing.assert_allclose(model.airfoil[:, 0, :],
                               np.repeat(angles[:, None], 2, axis=1))
    np.testing.assert_allclose(model.airfoil[:, 1, 0], angles / 20)
    np.testing.assert_allclose(model.airfoil[:, 2, 0], .2 + angles / 200)
    np.testing.assert_allclose(model.airfoil[:, 3, 0], angles / 100)
    np.testing.assert_allclose(model.airfoil[:, 1, 1], .5 + angles / 40)
    np.testing.assert_allclose(model.airfoil[:, 2, 1], .05 + angles / 400)

    # At a 10-degree inflow angle, the first midpoint has 3 degrees of
    # twist. Its interpolated polar therefore uses 7 degrees of attack.
    phi = np.deg2rad(10)
    lift, drag = 7 / 20, .2 + 7 / 200
    normal, tangent, moment = model._coefficients(phi, 0, 0)
    np.testing.assert_allclose(
        [normal, tangent, moment],
        [lift*np.cos(phi) + drag*np.sin(phi),
         lift*np.sin(phi) - drag*np.cos(phi), 7 / 100],
        rtol=0, atol=1e-14,
    )

    # The public source factory must use the same blade reader while placing
    # the published tilted shaft and hub in the platform body frame.
    rotor = MostRotor.from_iea15mw(tmp_path)
    np.testing.assert_array_equal(rotor.bem.radius, model.radius)
    tilt = np.deg2rad(6)
    np.testing.assert_allclose(rotor.shaft_axis,
                               [np.cos(tilt), 0, -np.sin(tilt)],
                               rtol=0, atol=1e-15)
    np.testing.assert_allclose(
        rotor.hub_offset,
        [-12.098 * np.cos(tilt), 0,
         15 + 129.386 + 4.3495 + 12.098 * np.sin(tilt) + 14.4],
        rtol=0, atol=1e-12,
    )
    assert rotor.inertia > 0


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


def test_short_coupled_platform_accounts_for_tower_mass_and_mooring_load():
    """Advance real BEM, controller, tower, mooring, and platform paths together."""
    origin = np.zeros(3)
    # Collocated masses make the initial heave balance independently
    # calculable while retaining all turbine reaction components.
    properties = {
        "tower": {"offset": 0, "mass": 100, "cog": origin,
                  "Inertia": origin},
        "nacelle": {"mass_yawBearing": 2, "cog_yawBearing": origin,
                    "mass": 3, "cog": origin, "Inertia": origin,
                    "tiltangle": 0},
        "hub": {"mass": 4, "cog": origin, "Rhub": 0, "precone": 0,
                "Inertia": origin},
        "blade": {"mass": 5, "cog_rel": origin, "Inertia": origin},
        "Inertia_Rotor_cogHub": 1,
    }
    tower = MostTowerReaction(properties, origin)
    platform = MostPlatformHydrodynamics(
        mass=1e8, inertia=np.eye(3) * 1e9,
        equilibrium_pose=np.zeros(6),
        restoring_stiffness=np.zeros((6, 6)),
        static_restoring=np.zeros(6), drag=np.zeros((6, 6)),
        radiation_time=np.array([0., .05, .1]),
        radiation_kernel=np.zeros((3, 6, 6)),
        added_mass=np.zeros((6, 6)), tower_base_height=0,
    )
    rotor = MostRotor(
        bem=_small_rotor(), controller=MostBaselineController.iea15mw(),
        inertia=1e9, hub_offset=origin,
        shaft_axis=np.array([1., 0., 0.]),
    )
    coupled = MostCoupled(platform, rotor, tower)
    time = np.array([0., .05, .1])
    waves = np.zeros((len(time), 6))
    wind = MostConstantWind(10, end_time=time[-1])
    converged = coupled.simulate(time, waves, wind)
    causal = coupled.simulate_causal(time, waves, wind)

    mooring_heave = MostStaticMooring().force(np.zeros(6))[0][2]
    turbine_mass = 100 + 2 + 3 + 4 + 3 * 5
    expected_initial_heave_acceleration = (
        mooring_heave - turbine_mass * tower.gravity
    ) / (platform.mass + turbine_mass)
    np.testing.assert_allclose(
        converged.platform.acceleration[0, 2],
        expected_initial_heave_acceleration, rtol=0, atol=1e-9,
    )
    assert converged.iterations >= 2
    assert converged.position_residual < 1e-6
    assert converged.velocity_residual < 1e-6
    assert converged.platform.position[-1, 0] > 0  # wind thrust
    assert converged.platform.position[-1, 2] < 0  # mooring and weight
    assert np.all(np.diff(converged.rotor.rotor_speed) < 0)
    np.testing.assert_allclose(
        causal.platform.position, converged.platform.position,
        rtol=0, atol=1e-7,
    )
    np.testing.assert_allclose(
        causal.rotor.rotor_speed, converged.rotor.rotor_speed,
        rtol=0, atol=1e-7,
    )

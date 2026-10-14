"""Check the floating OWC joint against a solvable two-mass limit."""

from types import SimpleNamespace

import numpy as np

from wecsim import MoorDyn, RegularWave, WEC, solve_floating_owc
import wecsim.dynamics.floatingOwcDynamics as floating_owc


class _ZeroMooring:
    def __init__(self):
        self.calls = 0
        self.closed = False

    def start(self, pose, velocity):
        np.testing.assert_array_equal(pose, np.zeros(6))
        np.testing.assert_array_equal(velocity, np.zeros(6))
        return self

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.closed = True

    def step(self, pose, velocity, time, dt):
        self.calls += 1
        assert np.isfinite(pose).all() and np.isfinite(velocity).all()
        assert dt > 0 and time >= 0
        return np.zeros(6)


class _NoAirForce:
    def force_on_column(self, pressure):
        return 0.0

    def pressure_derivative(self, pressure, displacement, speed, rotor_speed):
        return 0.0


class _SteadyTurbine:
    def evaluate(self, pressure, rotor_speed):
        zeros = np.zeros_like(np.asarray(rotor_speed, dtype=float))
        return SimpleNamespace(speed_derivative=zeros, load_power=zeros,
                               pneumatic_power=zeros)


def _two_mass_hydro(masses, floater_excitation=0):
    """A two-body limit with prescribed vertical excitation and no hydro memory."""
    def body_data(_h5, number, *_):
        mass = masses[number - 1]
        excitation = np.zeros(6)
        if number == 1:
            excitation[2] = floater_excitation
        return {
            "center": np.array([0.0, 0.0, float(number - 1)]),
            "mass": np.diag([mass] * 6),
            "damping": np.zeros((6, 6)),
            "restoring": np.zeros((6, 6)),
            "restoring_at_equilibrium": np.zeros(6),
            "vertical_bias": 0.0,
            "excitation_real": excitation,
            "excitation_imag": np.zeros(6),
        }
    return body_data


def test_coupled_slider_conserves_momentum_and_loses_only_pto_energy(
        tmp_path, monkeypatch):
    """The two body solver should reduce to an analytic damped slider."""
    floater_mass, column_mass, damping = 1000.0, 100.0, 50.0
    masses = (floater_mass, column_mass)

    monkeypatch.setattr(floating_owc, "_body_data", _two_mass_hydro(masses))
    library, lines = tmp_path / "libmoordyn.so", tmp_path / "lines.txt"
    library.touch()
    lines.touch()
    mooring = MoorDyn(library, lines)
    session = _ZeroMooring()
    monkeypatch.setattr(mooring, "start", session.start)
    monkeypatch.setattr(mooring, "step", session.step)

    initial = np.zeros(16)
    initial[13] = 0.2  # Relative column/floater heave speed.
    initial[15] = 150
    response = solve_floating_owc(
        library, mooring, dt=.01, end_time=1, wave_height=0,
        pto_damping=damping, chamber=_NoAirForce(), turbine=_SteadyTurbine(),
        initial_state=initial,
    )

    reduced_mass = floater_mass * column_mass / (floater_mass + column_mass)
    decay = damping / reduced_mass
    relative_speed = initial[13] * np.exp(-decay * response.time)
    stroke = initial[13] * (1 - np.exp(-decay * response.time)) / decay
    total_momentum = column_mass * initial[13]
    floater_speed = (total_momentum - column_mass * relative_speed) / (
        floater_mass + column_mass)

    assert session.calls == 100 and session.closed
    np.testing.assert_allclose(response.stroke, stroke, rtol=0, atol=1e-11)
    np.testing.assert_allclose(response.stroke_speed, relative_speed,
                               rtol=0, atol=1e-11)
    np.testing.assert_allclose(response.floater_velocity[:, 2], floater_speed,
                               rtol=0, atol=1e-11)
    np.testing.assert_allclose(response.column_velocity[:, 2],
                               floater_speed + relative_speed,
                               rtol=0, atol=1e-11)
    np.testing.assert_allclose(response.column_pose[:, 2],
                               response.floater_pose[:, 2] + 1 + stroke,
                               rtol=0, atol=1e-11)
    np.testing.assert_allclose(
        floater_mass * response.floater_velocity[:, 2]
        + column_mass * response.column_velocity[:, 2],
        total_momentum, rtol=0, atol=1e-11,
    )
    kinetic = (floater_mass * response.floater_velocity[:, 2]**2
               + column_mass * response.column_velocity[:, 2]**2) / 2
    expected_loss = reduced_mass * initial[13]**2 / 2 * (
        1 - np.exp(-2 * decay * response.time))
    np.testing.assert_allclose(kinetic[0] - kinetic, expected_loss,
                               rtol=0, atol=1e-11)
    np.testing.assert_allclose(response.pto_force,
                               -damping * response.stroke_speed,
                               rtol=0, atol=1e-12)
    np.testing.assert_allclose(response.pto_dissipated_power,
                               damping * response.stroke_speed**2,
                               rtol=0, atol=1e-12)


def test_public_floating_owc_routes_wave_motion_and_pto_work(tmp_path, monkeypatch):
    """The public builder must preserve slider geometry, force, and energy."""
    monkeypatch.setattr(
        floating_owc, "_body_data", _two_mass_hydro((1000.0, 100.0), 1000.0),
    )
    hydro = tmp_path / "owc.h5"
    library = tmp_path / "libmoordyn.so"
    lines = tmp_path / "lines.txt"
    for path in (hydro, library, lines):
        path.touch()
    mooring = MoorDyn(library, lines)
    session = _ZeroMooring()
    monkeypatch.setattr(mooring, "start", session.start)
    monkeypatch.setattr(mooring, "step", session.step)

    wec = WEC("Two-mass OWC through the public API")
    floater = wec.body("floater", hydro, hydro_body=1, mass=1000,
                       inertia=(1, 1, 1))
    column = wec.body("column", hydro, hydro_body=2, mass=100,
                      inertia=(1, 1, 1))
    wec.floating_owc(floater, column, moordyn=mooring,
                     column_height=1, column_diameter=1, pto_damping=50)
    result = wec.run(RegularWave(.5, 8), dt=.01, end_time=.5, ramp_time=0)

    assert result.time.shape == (51,)
    assert session.calls == 50 and session.closed
    np.testing.assert_allclose(
        result.wave_elevation, .25 * np.cos(2 * np.pi * result.time / 8),
        rtol=0, atol=1e-15,
    )
    floater_state = result.bodies["floater"]
    column_state = result.bodies["column"]
    pto = result.ptos["column_slider"]
    np.testing.assert_allclose(
        pto.stroke, column_state.position[:, 2] - floater_state.position[:, 2] - 1,
        rtol=0, atol=1e-12,
    )
    np.testing.assert_allclose(
        pto.velocity, column_state.velocity[:, 2] - floater_state.velocity[:, 2],
        rtol=0, atol=1e-12,
    )
    np.testing.assert_allclose(pto.force, -50 * pto.velocity, rtol=0, atol=1e-12)
    np.testing.assert_allclose(pto.absorbed_power, 50 * pto.velocity**2,
                               rtol=0, atol=1e-12)
    assert np.ptp(pto.stroke) > .02
    absorbed_work = np.sum(
        .5 * (pto.absorbed_power[1:] + pto.absorbed_power[:-1])
        * np.diff(result.time)
    )
    assert absorbed_work > .05
    extras = dict(result.raw.extra_outputs)
    assert np.isfinite(extras["chamber_pressure"]).all()
    assert np.isfinite(extras["turbine_speed"]).all()
    np.testing.assert_array_equal(extras["moordyn_connection_force"],
                                  np.zeros((51, 6)))

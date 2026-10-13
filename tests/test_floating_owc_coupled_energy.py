"""Check the floating OWC joint against a solvable two-mass limit."""

from types import SimpleNamespace

import numpy as np

from wecsim import MoorDyn, solve_floating_owc
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


def test_coupled_slider_conserves_momentum_and_loses_only_pto_energy(
        tmp_path, monkeypatch):
    """The two body solver should reduce to an analytic damped slider."""
    floater_mass, column_mass, damping = 1000.0, 100.0, 50.0
    masses = (floater_mass, column_mass)

    def body_data(_h5, number, *_):
        mass = masses[number - 1]
        return {
            "center": np.array([0.0, 0.0, float(number - 1)]),
            "mass": np.diag([mass] * 3 + [mass] * 3),
            "damping": np.zeros((6, 6)),
            "restoring": np.zeros((6, 6)),
            "restoring_at_equilibrium": np.zeros(6),
            "vertical_bias": 0.0,
            "excitation_real": np.zeros(6),
            "excitation_imag": np.zeros(6),
        }

    monkeypatch.setattr(floating_owc, "_body_data", body_data)
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

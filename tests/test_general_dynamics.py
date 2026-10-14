"""Check the shared dynamics stage against independent oscillator solutions."""

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from wecsim.dynamics.generalDynamics import (
    BodyMotion, DynamicBody, GeneralizedDynamics,
)
from wecsim.hydrodynamics.radiation import RadiationStateSpace


def _oscillator(*, mass=2.0, stiffness=8.0, radiation=None,
                pto_damping=0.0, pto_stiffness=0.0, pto_equilibrium=0.0,
                radiation_discretization="trapezoid", added_mass=0.0,
                added_mass_delay=None, nonlinear_force=None):
    jacobian = np.zeros((6, 1))
    jacobian[2, 0] = 1
    rigid_mass = np.zeros((6, 6))
    rigid_mass[2, 2] = mass
    restoring = np.zeros((6, 6))
    restoring[2, 2] = stiffness
    added = np.zeros((6, 6))
    added[2, 2] = added_mass

    def motion(q, v):
        displacement = np.zeros(6)
        displacement[2] = q[0]
        return BodyMotion(displacement, jacobian, np.zeros(6))

    body = DynamicBody(
        rigid_mass=rigid_mass,
        added_mass=(added,),
        damping=(np.zeros((6, 6)),),
        restoring=restoring,
        static_force=np.zeros(6),
        reference_position=np.zeros(6),
        motion=motion,
        excitation=lambda t: np.zeros(6),
        radiation_kernel=radiation,
    )
    return GeneralizedDynamics(
        (body,), 1, pto_damping=np.array([[pto_damping]]),
        pto_stiffness=np.array([[pto_stiffness]]),
        pto_equilibrium=np.array([pto_equilibrium]),
        radiation_discretization=radiation_discretization,
        added_mass_delay=added_mass_delay,
        nonlinear_force=nonlinear_force,
    )


def test_regular_rk4_matches_damped_oscillator_solution():
    mass, stiffness, damping = 2.0, 8.0, 0.5
    system = _oscillator(mass=mass, stiffness=stiffness, pto_damping=damping)
    response = system.integrate(
        dt=0.01, end_time=10, initial_coordinate=np.array([1.0]),
    )
    decay = damping / (2 * mass)
    frequency = np.sqrt(stiffness / mass - decay**2)
    expected = np.exp(-decay * response.time) * (
        np.cos(frequency * response.time)
        + decay / frequency * np.sin(frequency * response.time)
    )
    assert np.max(np.abs(response.coordinate[:, 0] - expected)) < 1e-7
    assert np.max(np.abs(response.body_position[:, 0, 2] - expected)) < 1e-7


def test_pto_equilibrium_shift_changes_force_and_motion():
    system = _oscillator(stiffness=0, pto_stiffness=8,
                         pto_equilibrium=1)
    response = system.integrate(dt=0.01, end_time=2)
    expected = 1 - np.cos(2 * response.time)
    assert np.max(np.abs(response.coordinate[:, 0] - expected)) < 1e-8
    np.testing.assert_allclose(
        system.acceleration(0, np.array([0.0]), np.array([0.0])),
        [4.0], rtol=0, atol=1e-12,
    )


def test_radiation_memory_matches_augmented_state_equation():
    dt, duration, mass, stiffness = 0.01, 5.0, 2.0, 8.0
    strength, decay = 1.5, 0.8
    lag = np.arange(round(duration / dt) + 1) * dt
    kernel = np.zeros((len(lag), 6, 6))
    kernel[:, 2, 2] = strength * np.exp(-decay * lag)
    system = _oscillator(mass=mass, stiffness=stiffness, radiation=kernel)
    response = system.integrate(
        dt=dt, end_time=duration, initial_coordinate=np.array([1.0]),
    )

    def derivative(t, state):
        q, v, radiation = state
        return [v, (-stiffness * q - radiation) / mass,
                strength * v - decay * radiation]

    expected = solve_ivp(
        derivative, (0, duration), [1, 0, 0], t_eval=response.time,
        rtol=1e-11, atol=1e-13,
    )
    assert expected.success
    assert np.max(np.abs(response.coordinate[:, 0] - expected.y[0])) < 3e-4
    assert np.max(np.abs(response.speed[:, 0] - expected.y[1])) < 6e-4


def test_two_body_fitted_radiation_matches_independent_state_equation():
    """The fitted force on body two must project through its reversed motion."""
    mass, stiffness, strength, decay = 5.0, 10.0, 1.5, 0.8

    def body(sign, rigid_mass):
        jacobian = np.zeros((6, 1))
        jacobian[2, 0] = sign
        inertia = np.zeros((6, 6))
        inertia[2, 2] = rigid_mass
        restoring = np.zeros((6, 6))
        restoring[2, 2] = stiffness / 2

        def motion(q, v):
            return BodyMotion(jacobian[:, 0] * q[0], jacobian, np.zeros(6))

        return DynamicBody(
            rigid_mass=inertia,
            added_mass=(np.zeros((6, 6)), np.zeros((6, 6))),
            damping=(np.zeros((6, 6)), np.zeros((6, 6))),
            restoring=restoring,
            static_force=np.zeros(6),
            reference_position=np.zeros(6),
            motion=motion,
            excitation=lambda t: np.zeros(6),
        )

    # The state is driven by body one's +q velocity. Body two moves as -q;
    # a negative physical force on body two becomes positive generalized
    # radiation, which the dynamics subtracts from the oscillator balance.
    b = np.zeros((1, 12))
    b[0, 2] = strength
    c = np.zeros((12, 1))
    c[8, 0] = -1
    radiation = RadiationStateSpace(np.array([[-decay]]), b, c)
    system = GeneralizedDynamics((body(1, 2), body(-1, 3)), 1)
    response = system.integrate(
        dt=0.01, end_time=5, initial_coordinate=np.array([1.0]),
        radiation_state_space=radiation,
    )

    expected = solve_ivp(
        lambda t, state: [
            state[1], (-stiffness * state[0] - state[2]) / mass,
            strength * state[1] - decay * state[2],
        ],
        (0, 5), [1, 0, 0], t_eval=response.time,
        rtol=1e-11, atol=1e-13,
    )
    assert expected.success
    np.testing.assert_allclose(response.coordinate[:, 0], expected.y[0],
                               rtol=0, atol=1e-7)
    np.testing.assert_allclose(response.speed[:, 0], expected.y[1],
                               rtol=0, atol=1e-7)
    np.testing.assert_allclose(response.acceleration[:, 0],
                               (-stiffness * expected.y[0] - expected.y[2]) / mass,
                               rtol=0, atol=1e-7)
    np.testing.assert_allclose(response.body_position[:, 0, 2], expected.y[0],
                               rtol=0, atol=1e-7)
    np.testing.assert_allclose(response.body_position[:, 1, 2], -expected.y[0],
                               rtol=0, atol=1e-7)


def test_discontinuous_sampled_force_uses_accepted_selection_only_on_failure():
    """An implicit endpoint can have no solution across a force jump."""
    from dataclasses import replace

    class JumpForce:
        def __init__(self):
            self.fallback_times = []

        def __call__(self, at_time, coordinate, speed):
            force = np.zeros(6)
            force[2] = 1 if coordinate[0] < .25 else -1
            return force

        def freeze_step_excitation(self, at_time, accepted_coordinate,
                                   trial_coordinate, next_coordinate):
            if ((trial_coordinate[0] < .25)
                    == (next_coordinate[0] < .25)):
                return None
            self.fallback_times.append(at_time)
            return self(at_time, accepted_coordinate, np.zeros(1))

    kernel = np.zeros((2, 6, 6))
    system = _oscillator(mass=1, stiffness=0, radiation=kernel)
    source = JumpForce()
    system.bodies = (replace(system.bodies[0], state_excitation=source),)
    response = system.integrate(dt=1, end_time=1)
    np.testing.assert_allclose(response.coordinate[:, 0], [0, .5])
    np.testing.assert_allclose(response.speed[:, 0], [0, 1])
    np.testing.assert_allclose(response.acceleration[:, 0], [1, -1])
    assert source.fallback_times == [1]

    no_fallback = _oscillator(mass=1, stiffness=0, radiation=kernel)
    no_fallback.bodies = (replace(no_fallback.bodies[0],
                                  state_excitation=lambda t, q, v: source(t, q, v)),)
    with pytest.raises(RuntimeError, match="did not converge"):
        no_fallback.integrate(dt=1, end_time=1)


def test_trapezoidal_radiation_halves_full_window_endpoint():
    kernel = np.zeros((3, 6, 6))
    kernel[-1, 2, 2] = 1
    system = _oscillator(mass=1, stiffness=0, radiation=kernel)
    response = system.integrate(
        dt=0.1, end_time=0.2, initial_speed=np.array([1.0]),
    )
    # At t=0.2 s the oldest velocity first reaches the 0.2 s memory
    # boundary. The trapezoidal convolution gives that endpoint half weight.
    np.testing.assert_allclose(response.acceleration[:, 0],
                               [0, 0, -0.05], rtol=0, atol=1e-12)


def test_stateful_external_force_advances_once_per_memory_step():
    kernel = np.zeros((3, 6, 6))
    system = _oscillator(mass=2, stiffness=0, radiation=kernel)
    calls = []

    def step_force(time, dt, predicted_coordinate, predicted_speed):
        calls.append((time, dt, predicted_coordinate.copy(),
                      predicted_speed.copy()))
        return np.array([2.0])

    response = system.integrate(dt=0.1, end_time=0.2,
                                step_force=step_force)
    assert len(calls) == 2
    np.testing.assert_allclose([call[0] for call in calls], [0, 0.1])
    np.testing.assert_allclose([call[1] for call in calls], [0.1, 0.1])
    np.testing.assert_allclose(calls[0][2:], [[0], [0]])
    np.testing.assert_allclose(response.acceleration[:, 0], [0, 1, 1],
                               rtol=0, atol=1e-12)
    np.testing.assert_allclose(response.speed[:, 0], [0, 0.05, 0.15],
                               rtol=0, atol=1e-12)
    np.testing.assert_allclose(response.coordinate[:, 0],
                               [0, 0.0025, 0.0125], rtol=0, atol=1e-12)


def test_delayed_added_mass_matches_independent_oscillator_recurrence():
    kernel = np.zeros((3, 6, 6))
    dt, mass, added, stiffness = 0.1, 2.0, 1.0, 8.0
    delay = 1e-7
    system = _oscillator(
        mass=mass, stiffness=stiffness, radiation=kernel,
        added_mass=added, added_mass_delay=delay,
    )
    response = system.integrate(
        dt=dt, end_time=0.3, initial_coordinate=np.array([1.0]),
    )
    q = np.zeros(4)
    v = np.zeros(4)
    a = np.zeros(4)
    q[0] = 1
    adjusted_mass = mass + 2 * added
    applied_mass = added - 2 * added
    a[0] = -stiffness * q[0] / adjusted_mass
    for step in range(1, 4):
        delayed = (a[0] if step == 1 else
                   a[step - 1] + (1 - delay / dt)
                   * (a[step - 1] - a[step - 2]))
        a[step] = (
            -stiffness * (q[step - 1] + dt * v[step - 1]
                          + dt**2 / 4 * a[step - 1])
            - applied_mass * delayed
        ) / (adjusted_mass + stiffness * dt**2 / 4)
        v[step] = v[step - 1] + dt * (a[step - 1] + a[step]) / 2
        q[step] = q[step - 1] + dt * (v[step - 1] + v[step]) / 2
    np.testing.assert_allclose(response.coordinate[:, 0], q, rtol=0, atol=1e-12)
    np.testing.assert_allclose(response.speed[:, 0], v, rtol=0, atol=1e-12)
    np.testing.assert_allclose(response.acceleration[:, 0], a, rtol=0, atol=1e-12)
    repeated = system.integrate(
        dt=dt, end_time=0.3, initial_coordinate=np.array([1.0]),
    )
    np.testing.assert_array_equal(repeated.coordinate, response.coordinate)


def test_implicit_added_mass_preserves_undamped_oscillator_energy():
    kernel = np.zeros((101, 6, 6))
    mass, added, stiffness = 2.0, 1.0, 8.0
    system = _oscillator(mass=mass, stiffness=stiffness, radiation=kernel,
                         added_mass=added)
    response = system.integrate(
        dt=0.01, end_time=1, initial_coordinate=np.array([1.0]),
    )
    np.testing.assert_allclose(
        response.acceleration[0, 0], -stiffness / (mass + added),
        rtol=0, atol=1e-12,
    )
    energy = (0.5 * (mass + added) * response.speed[:, 0]**2
              + 0.5 * stiffness * response.coordinate[:, 0]**2)
    assert np.max(np.abs(energy - energy[0])) < 1e-10


def test_adaptive_unilateral_stop_conserves_energy_and_rebounds():
    mass, stiffness, bound = 2.0, 10_000.0, 0.1

    def stop(coordinate, speed):
        del speed
        return np.array([-stiffness * max(coordinate[0] - bound, 0)])

    system = _oscillator(mass=mass, stiffness=0,
                         nonlinear_force=stop)
    response = system.integrate(
        dt=0.02, end_time=0.8, initial_speed=np.array([1.0]),
        adaptive_regular=True,
    )
    position, speed = response.coordinate[:, 0], response.speed[:, 0]
    expected_peak = bound + np.sqrt(mass / stiffness)
    assert abs(np.max(position) - expected_peak) < 2e-4
    assert speed[-1] < -0.999
    energy = (0.5 * mass * speed**2
              + 0.5 * stiffness * np.maximum(position - bound, 0)**2)
    assert np.max(np.abs(energy - energy[0])) < 2e-5


def test_fir_uses_full_sampled_taps_and_holds_force_during_each_step():
    kernel = np.zeros((2, 6, 6))
    kernel[:, 2, 2] = [2, 3]
    system = _oscillator(
        mass=1, stiffness=0, radiation=kernel,
        radiation_discretization="fir",
    )
    response = system.integrate(
        dt=0.1, end_time=0.2, initial_speed=np.array([1.0]),
    )
    np.testing.assert_allclose(response.speed[:, 0],
                               [1, 0.98, 0.9304], rtol=0, atol=1e-12)
    np.testing.assert_allclose(response.coordinate[:, 0],
                               [0, 0.099, 0.19452], rtol=0, atol=1e-12)
    np.testing.assert_allclose(response.acceleration[:, 0],
                               [-0.2, -0.496, -0.48008],
                               rtol=0, atol=1e-12)

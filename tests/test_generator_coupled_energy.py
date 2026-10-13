"""Energy and momentum guard for a generator coupled to two free bodies."""

import numpy as np
from scipy.integrate import trapezoid

from wecsim.components.directLinearGenerator import DirectLinearGenerator
from wecsim.dynamics.generalDynamics import (
    BodyMotion, DynamicBody, GeneralizedDynamics,
)
from wecsim.dynamics.ptoConnections import LinearPTOConnection


def _free_heave_body(index, mass):
    jacobian = np.zeros((6, 2))
    jacobian[2, index] = 1

    def motion(coordinate, _speed):
        displacement = np.zeros(6)
        displacement[2] = coordinate[index]
        return BodyMotion(displacement, jacobian, np.zeros(6))

    rigid_mass = np.eye(6)
    rigid_mass[2, 2] = mass
    zero_matrix = np.zeros((6, 6))
    return DynamicBody(
        rigid_mass=rigid_mass,
        added_mass=(zero_matrix, zero_matrix),
        damping=(zero_matrix, zero_matrix),
        restoring=zero_matrix,
        static_force=np.zeros(6),
        reference_position=np.zeros(6),
        motion=motion,
        excitation=lambda _time: np.zeros(6),
    )


def test_two_body_generator_conserves_momentum_and_accounts_for_energy():
    """Internal PTO force transfers momentum and dissipates measured energy."""
    generator = DirectLinearGenerator(
        stator_resistance=4.58, friction=-100,
        pole_pitch=.072, magnet_flux=8, inductance=.285,
        load_resistance=-117.6471,
    )
    connection = LinearPTOConnection(
        "generator", np.array([1., -1.]), damping=0,
        stiffness=0, equilibrium_position=0,
    )
    dynamics = GeneralizedDynamics(
        (_free_heave_body(0, 2), _free_heave_body(1, 3)), 2,
        linear_generators=((connection, generator),),
    )
    result = dynamics.integrate(
        dt=.0001, end_time=.05, initial_speed=np.array([.4, 0.]),
    )
    assert result.linear_generator_state.shape == (501, 1, 3)
    assert result.linear_generator_force.shape == (501, 1)
    np.testing.assert_allclose(result.linear_generator_force[0, 0], -40,
                               rtol=0, atol=1e-12)
    np.testing.assert_allclose(result.acceleration[0], [-20, 40 / 3],
                               rtol=0, atol=1e-12)
    momentum = 2 * result.speed[:, 0] + 3 * result.speed[:, 1]
    np.testing.assert_allclose(momentum, .8, rtol=0, atol=1e-12)

    # Mechanical plus magnetic energy can only leave through friction,
    # stator resistance, and the passive external electrical load.
    relative_speed = result.speed[:, 0] - result.speed[:, 1]
    flux = result.linear_generator_state[:, 0, :]
    current_squared = (
        (flux[:, 0] - generator.magnet_flux)**2 + flux[:, 1]**2
    ) / generator.inductance**2
    kinetic = (2 * result.speed[:, 0]**2 + 3 * result.speed[:, 1]**2) / 2
    magnetic = (
        (flux[:, 0] - generator.magnet_flux)**2 + flux[:, 1]**2
    ) / (2 * generator.inductance)
    friction_loss = -generator.friction * relative_speed**2
    electrical_loss = (
        generator.stator_resistance - generator.load_resistance
    ) * current_squared
    assert np.all(friction_loss >= 0) and np.all(electrical_loss >= 0)
    lost = kinetic[0] + magnetic[0] - kinetic[-1] - magnetic[-1]
    dissipated = trapezoid(friction_loss + electrical_loss, x=result.time)
    assert lost > 0
    np.testing.assert_allclose(lost, dissipated, rtol=0, atol=6e-6)

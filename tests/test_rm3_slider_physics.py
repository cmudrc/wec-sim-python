"""Independent geometry and energy checks for the RM3 floating joint."""

import numpy as np
import pytest

from wecsim.dynamics.generalDynamics import DynamicBody, GeneralizedDynamics
from wecsim.reference_models.rm3Regular import _rm3_slider_motion


@pytest.mark.parametrize("body_index,lever", [(0, 1.7), (1, -2.3)])
def test_slider_jacobian_and_curvature_are_position_derivatives(body_index, lever):
    coordinate = np.array([0.4, -0.2, 0.3, 0.37])
    speed = np.array([0.15, -0.4, 0.2, -0.31])
    step = 1e-6
    motion = _rm3_slider_motion(
        coordinate, speed, body_index=body_index, lever=lever,
    )

    # Check the six-dimensional body map without using the solver's formulas
    # for its derivatives.
    columns = []
    for axis in range(4):
        offset = np.eye(4)[axis] * step
        plus = _rm3_slider_motion(
            coordinate + offset, speed, body_index=body_index, lever=lever,
        ).displacement
        minus = _rm3_slider_motion(
            coordinate - offset, speed, body_index=body_index, lever=lever,
        ).displacement
        columns.append((plus - minus) / (2 * step))
    np.testing.assert_allclose(motion.jacobian, np.column_stack(columns),
                               rtol=0, atol=1e-9)

    plus = _rm3_slider_motion(
        coordinate + step * speed, speed,
        body_index=body_index, lever=lever,
    ).jacobian
    minus = _rm3_slider_motion(
        coordinate - step * speed, speed,
        body_index=body_index, lever=lever,
    ).jacobian
    numerical_curvature = ((plus - minus) / (2 * step)) @ speed
    np.testing.assert_allclose(motion.bias_acceleration, numerical_curvature,
                               rtol=0, atol=1e-9)


def test_unforced_pitched_sliders_conserve_kinetic_energy():
    """A constrained joint must not create power through its geometry."""
    own_added = np.diag([0.4, 0, 0.5, 0, 0.1, 0])
    cross_added = np.diag([0.1, 0, 0.05, 0, 0.02, 0])
    masses = (np.diag([3.0, 3.0, 3.0, 0, 2.0, 0]),
              np.diag([4.0, 4.0, 4.0, 0, 3.0, 0]))
    levers = (1.7, -2.3)
    bodies = []
    for index in (0, 1):
        blocks = ((own_added, cross_added) if index == 0
                  else (cross_added, own_added))
        bodies.append(DynamicBody(
            rigid_mass=masses[index], added_mass=blocks,
            damping=(np.zeros((6, 6)), np.zeros((6, 6))),
            restoring=np.zeros((6, 6)), static_force=np.zeros(6),
            reference_position=np.zeros(6),
            motion=lambda q, v, index=index: _rm3_slider_motion(
                q, v, body_index=index, lever=levers[index],
            ),
            excitation=lambda time: np.zeros(6),
        ))
    response = GeneralizedDynamics(tuple(bodies), 4).integrate(
        dt=0.01, end_time=10,
        initial_coordinate=np.array([0.2, 0.3, -0.2, 0.4]),
        initial_speed=np.array([0.15, 0.35, -0.25, 0.2]),
    )
    velocity = response.body_velocity
    energy = np.zeros(len(response.time))
    for output in (0, 1):
        for source in (0, 1):
            block = (masses[output] + own_added if output == source
                     else cross_added)
            energy += 0.5 * np.einsum(
                "ti,ij,tj->t", velocity[:, output], block,
                velocity[:, source],
            )
    assert np.max(np.abs(energy - energy[0])) / energy[0] < 1e-8

"""Physical and input-boundary checks for named WEC coordinates."""

import numpy as np
import pytest

from wecsim.dynamics.linearCoordinates import (
    build_coordinate_maps,
    initial_coordinate,
)


CENTERS = [np.array([1.0, 2.0, 3.0]), np.array([0.0, -3.0, 1.0])]
BODIES = [{}, {}]
NAMES = ["float", "spar"]


def _maps(*definitions):
    return build_coordinate_maps(
        {"coordinates": list(definitions)}, BODIES, NAMES, CENTERS,
    )


def test_shared_roll_about_two_attachment_points_preserves_virtual_work():
    # The float's pivot is one metre above its center in body coordinates;
    # the spar's pivot is the world origin. Opposite scales couple their rolls.
    maps, names = _maps({
        "name": "relative_roll",
        "motions": [
            {"body": "float", "dof": "roll", "scale": 2,
             "pivot": {"point": [0, 0, 1]}},
            {"body": "spar", "dof": "roll", "scale": -1,
             "pivot": {"world": [0, 0, 0]}},
        ],
    })
    assert names == ("relative_roll",)
    np.testing.assert_allclose(maps[0][:, 0], [0, 2, 0, 2, 0, 0])
    np.testing.assert_allclose(maps[1][:, 0], [0, 1, 3, -1, 0, 0])

    # A generalized speed must do the same instantaneous work as the
    # resulting six-component body speeds under arbitrary applied loads.
    speed = 0.7
    loads = [np.array([4, -5, 6, 7, -8, 9]),
             np.array([-2, 3, -4, 5, -6, 7])]
    body_power = sum(np.dot(mapping[:, 0] * speed, load)
                     for mapping, load in zip(maps, loads))
    generalized_force = sum(np.dot(mapping[:, 0], load)
                            for mapping, load in zip(maps, loads))
    np.testing.assert_allclose(body_power, speed * generalized_force,
                               rtol=0, atol=1e-12)


def test_sparse_named_initial_state_defaults_other_coordinates_to_zero():
    names = ("float_heave", "spar_heave", "relative_roll")
    np.testing.assert_array_equal(
        initial_coordinate({"relative_roll": 0.2}, names, "initial state"),
        [0, 0, 0.2],
    )
    with pytest.raises(ValueError, match="unknown coordinates"):
        initial_coordinate({"pitch": 0.2}, names, "initial state")
    with pytest.raises(ValueError, match="one finite value"):
        initial_coordinate([0, np.nan, 0], names, "initial state")


@pytest.mark.parametrize("motion, message", [
    ({"body": True, "dof": "heave"}, "body is unknown"),
    ({"body": "float", "dof": "heave", "scale": 0}, "scale must be finite"),
    ({"body": "float", "dof": "heave", "scale": np.inf},
     "scale must be finite"),
    ({"body": "float", "dof": "heave", "pivot": {"world": [0, 0, 0]}},
     "only rotational motions"),
    ({"body": "float", "dof": "roll", "pivot": {"point": [0, 0]}},
     "finite three-vector"),
    ({"body": "float", "dof": "roll", "pivot": {"body": [0, 0, 0]}},
     "world or body-local"),
])
def test_unsafe_coordinate_motion_is_rejected(motion, message):
    with pytest.raises(ValueError, match=message):
        _maps({"name": "motion", "motions": [motion]})


def test_duplicate_body_motion_cannot_double_count_kinematics():
    with pytest.raises(ValueError, match="cannot repeat"):
        _maps({"name": "heave", "motions": [
            {"body": "float", "dof": "heave"},
            {"body": 1, "dof": "heave"},
        ]})


def test_named_coordinates_cannot_mix_with_raw_body_maps():
    with pytest.raises(ValueError, match="named coordinates or"):
        build_coordinate_maps(
            {"coordinates": [{"name": "heave", "motions": [
                {"body": "float", "dof": "heave"}]}]},
            [{"coordinate_map": [[1], [0], [0], [0], [0], [0]]}, {}],
            NAMES, CENTERS,
        )

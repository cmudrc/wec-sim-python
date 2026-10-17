"""Keep Boolean flags out of physical initial displacement and speed."""

from pathlib import Path

import numpy as np
import pytest

from wecsim import NoWave, WEC
from wecsim.dynamics.linearCoordinates import initial_coordinate


SPHERE = Path(__file__).resolve().parents[1] / "examples/data/sphere.h5"


@pytest.mark.parametrize(("field", "value"), [
    ("initial_coordinate", {"heave": True}),
    ("initial_coordinate", [False]),
    ("initial_speed", {"heave": False}),
    ("initial_speed", [True]),
])
def test_public_wec_rejects_boolean_initial_states(field, value):
    wec = WEC("Sphere initial-state validation")
    sphere = wec.body("sphere", SPHERE)
    wec.coordinate("heave", sphere.move("heave"))
    with pytest.raises(ValueError, match="one finite value per coordinate"):
        wec.run(NoWave(), dt=.1, end_time=.2, radiation_memory=.2,
                **{field: value})


def test_numpy_boolean_initial_state_is_also_rejected():
    with pytest.raises(ValueError, match="one finite value per coordinate"):
        initial_coordinate(np.array([True]), ("heave",), "initial state")

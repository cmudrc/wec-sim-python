"""A body drag setting must stay attached to its body in a coupled WEC."""

from pathlib import Path

import numpy as np

from wecsim import NoWave, WEC


RM3 = Path(__file__).resolve().parents[1] / "examples/data/rm3.h5"


def _run(*, drag_form):
    wec = WEC("RM3 independent heaves", body_to_body=True)
    if drag_form == "vectors":
        drag = dict(quadratic_drag_coefficient=(0, 0, 1, 0, 0, 0),
                    quadratic_drag_area=(0, 0, 1000, 0, 0, 0))
    elif drag_form == "matrix":
        # D = rho * Cd * area / 2 at the default 1000 kg/m^3 density.
        drag = dict(quadratic_drag=np.diag((0, 0, 500_000, 0, 0, 0)))
    else:
        drag = {}
    float_body = wec.body("float", RM3, **drag)
    spar = wec.body("spar", RM3)
    wec.coordinate("float_heave", float_body.move("heave"))
    wec.coordinate("spar_heave", spar.move("heave"))
    return wec.run(
        NoWave(), dt=.01, end_time=.2, radiation_memory=.2,
        initial_speed={"float_heave": .2, "spar_heave": -.1},
    )


def test_one_body_drag_vectors_match_matrix_with_cross_body_hydrodynamics():
    vectors = _run(drag_form="vectors")
    matrix = _run(drag_form="matrix")
    control = _run(drag_form="none")

    for name in ("float", "spar"):
        for field in ("position", "velocity"):
            np.testing.assert_array_equal(
                getattr(vectors.bodies[name], field),
                getattr(matrix.bodies[name], field),
            )
    assert np.max(np.abs(vectors.bodies["float"].velocity
                         - control.bodies["float"].velocity)) > 1e-3
    # The spar changes through BEM body-to-body coupling even though it has
    # no drag of its own. This checks the assembled force, not just the input.
    assert np.max(np.abs(vectors.bodies["spar"].velocity
                         - control.bodies["spar"].velocity)) > 1e-5

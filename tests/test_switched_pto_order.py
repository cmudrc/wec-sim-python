"""A switched PTO and a passive PTO must not depend on declaration order."""

from pathlib import Path

import numpy as np

from wecsim import DeclutchingControl, RegularWave, WEC, WorldPoint


SPHERE = Path(__file__).resolve().parents[1] / "examples/data/sphere.h5"


def _run(order):
    wec = WEC("Sphere with two PTOs")
    sphere = wec.body("sphere", SPHERE)
    wec.coordinate("heave", sphere.move("heave"))
    for name in order:
        if name == "switched":
            wec.pto(
                name, WorldPoint(0, 0, -5), sphere.at(0, 0, 0),
                control=DeclutchingControl(150_000, .12),
            )
        else:
            wec.pto(
                name, WorldPoint(0, 0, -5), sphere.at(0, 0, 0),
                stiffness=40_000, damping=30_000,
            )
    return wec.run(
        RegularWave(0, 8), dt=.01, end_time=10, ramp_time=0,
        initial_coordinate={"heave": .5},
    )


def test_switched_and_passive_pto_order_preserves_motion_and_work():
    switched_first = _run(("switched", "passive"))
    passive_first = _run(("passive", "switched"))

    for field in ("position", "velocity"):
        np.testing.assert_array_equal(
            getattr(switched_first.bodies["sphere"], field),
            getattr(passive_first.bodies["sphere"], field),
        )
    for name in ("switched", "passive"):
        for field in ("stroke", "velocity", "force", "absorbed_power"):
            np.testing.assert_array_equal(
                getattr(switched_first.ptos[name], field),
                getattr(passive_first.ptos[name], field),
            )

    switched = switched_first.ptos["switched"]
    passive = switched_first.ptos["passive"]
    # The free decay crosses several velocity reversals, exercising both
    # engaged damping and the controller's zero-force declutched interval.
    assert np.count_nonzero(switched.force == 0) > 10
    assert np.count_nonzero(switched.force != 0) > 10
    assert np.ptp(switched.stroke) > .5
    assert np.min(switched.absorbed_power) >= 0
    assert np.min(passive.absorbed_power) >= 0
    np.testing.assert_allclose(
        switched.absorbed_power, -switched.force * switched.velocity,
        rtol=0, atol=1e-9,
    )
    np.testing.assert_allclose(
        passive.absorbed_power, 30_000 * passive.velocity**2,
        rtol=0, atol=1e-9,
    )

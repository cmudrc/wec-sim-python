"""Reject misleading rigid-body inputs before a public WEC is mutated."""

from pathlib import Path

import numpy as np
import pytest

from wecsim import NoWave, WEC


SPHERE = Path(__file__).resolve().parents[1] / "examples/data/sphere.h5"


def test_moving_body_rejects_invalid_mass_and_inertia_then_runs():
    wec = WEC("rigid input recovery")
    for mass in (True, np.bool_(True), 0, -1, np.nan, np.inf, "250000", [250000]):
        with pytest.raises(ValueError, match="body mass.*positive"):
            wec.body("sphere", SPHERE, mass=mass)
        assert not wec.bodies
    for inertia in ((0, True, 0), (0, np.bool_(True), 0),
                    (0, -1, 0), (0, np.nan, 0), (0, 0), None):
        with pytest.raises(ValueError, match="body inertia.*nonnegative"):
            wec.body("sphere", SPHERE, inertia=inertia)
        assert not wec.bodies

    body = wec.body("sphere", SPHERE, mass=np.float64(250_000),
                    inertia=(0, np.int64(10), 0))
    assert body.mass == 250_000.0 and body.inertia == (0.0, 10.0, 0.0)
    wec.coordinate("heave", body.move("heave"))
    result = wec.run(NoWave(), dt=.05, end_time=.2,
                     initial_coordinate={"heave": .2})
    assert np.isfinite(result.bodies["sphere"].position).all()
    assert np.ptp(result.bodies["sphere"].position[:, 2]) > 0


def test_fixed_body_can_be_massless_but_rejects_invalid_rigid_inputs():
    wec = WEC("fixed input recovery")
    for settings in ({"mass": True}, {"mass": np.nan}, {"mass": -1},
                     {"inertia": (0, np.bool_(True), 0)},
                     {"inertia": (0, 0, np.inf)}):
        with pytest.raises(ValueError, match="body mass|body inertia"):
            wec.fixed_body("base", center_gravity=(0, 0, -10), **settings)
        assert not wec.bodies

    base = wec.fixed_body("base", center_gravity=(0, 0, -10),
                          mass=0, inertia=(0, 0, 0))
    assert base.mass == 0.0 and base.fixed is True

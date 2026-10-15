"""Rejected physical layouts must leave a Python WEC builder usable."""

from pathlib import Path

import numpy as np
import pytest

from wecsim import NoWave, RegularWave, WEC, WorldPoint


ROOT = Path(__file__).resolve().parents[1]
RM3 = ROOT / "examples/data/rm3.h5"
SPHERE = ROOT / "examples/data/sphere.h5"
OSWEC = ROOT / "examples/data/oswec.h5"


def test_rm3_rejected_unpaired_mooring_point_does_not_select_joint():
    wec = WEC("RM3")
    float_body = wec.body("float", RM3, inertia=(0, 21_306_090.66, 0))
    spar = wec.body("spar", RM3, inertia=(0, 94_407_091.24, 0))
    with pytest.raises(ValueError, match="needs a MoorDyn session"):
        wec.floating_joint(float_body, spar, moordyn_point=spar.at(0, 0, 21.5),
                           damping=1_200_000)

    wec.floating_joint(float_body, spar, damping=1_200_000)
    result = wec.run(RegularWave(.1, 8), dt=.1, end_time=.2, ramp_time=0)
    assert result.case["constraint"]["kind"] == "floating_joint"
    assert np.isfinite(result.bodies["float"].position).all()
    assert np.isfinite(result.ptos["relative_heave"].force).all()


def test_oswec_rejected_off_axis_hinge_does_not_select_layout():
    wec = WEC("OSWEC")
    flap = wec.body("flap", OSWEC, mass=127_000,
                    inertia=(0, 1.85e6, 0))
    base = wec.fixed_body("base", center_gravity=(0, 0, -10.9), mass=999)
    with pytest.raises(ValueError, match="z axis"):
        wec.fixed_hinge(flap, base, location=WorldPoint(1, 0, -10),
                        damping=12_000)

    wec.fixed_hinge(flap, base, damping=12_000)
    result = wec.run(RegularWave(.1, 8), dt=.1, end_time=.2, ramp_time=0)
    assert result.case["constraint"]["kind"] == "fixed_hinge"
    assert np.isfinite(result.bodies["flap"].position).all()
    assert np.isfinite(result.ptos["hinge"].force).all()


def test_sphere_rejected_foreign_pto_point_does_not_reserve_name():
    wec = WEC("Sphere")
    body = wec.body("sphere", SPHERE)
    wec.coordinate("heave", body.move("heave"))
    foreign = WEC("Other").body("sphere", SPHERE)
    with pytest.raises(ValueError, match="must belong to this WEC"):
        wec.pto("main", WorldPoint(0, 0, -5), foreign.at(0, 0, 0),
                damping=1_000)

    wec.pto("main", WorldPoint(0, 0, -5), body.at(0, 0, 0),
            damping=1_000)
    result = wec.run(NoWave(), dt=.01, end_time=.2,
                     radiation_memory=.2, initial_coordinate={"heave": .1})
    assert np.ptp(result.ptos["main"].stroke) > 0
    assert np.all(result.ptos["main"].absorbed_power >= 0)

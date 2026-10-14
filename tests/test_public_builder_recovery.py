"""Invalid public configurations must leave the WEC builder usable."""

from pathlib import Path

import pytest

from wecsim import (
    DeclutchingControl,
    DirectLinearGenerator,
    NoWave,
    SimpleDirectDrive,
    WEC,
    WorldPoint,
)


HYDRO = Path(__file__).resolve().parents[1] / "examples/data/rm3.h5"


def _heave_wec():
    wec = WEC("configured heave")
    body = wec.body("float", HYDRO)
    coordinate = wec.coordinate("heave", body.move("heave"))
    return wec, body, coordinate


def test_rejected_body_and_coordinate_do_not_corrupt_builder():
    wec, body, coordinate = _heave_wec()
    other = WEC("other")
    foreign = other.body("spar", HYDRO)

    with pytest.raises(ValueError, match="body name already exists"):
        wec.body("float", HYDRO)
    with pytest.raises(ValueError, match="coordinate motions must use bodies"):
        wec.coordinate("foreign", foreign.move("heave"))
    with pytest.raises(ValueError, match="coordinate name already exists"):
        wec.coordinate("heave", body.move("surge"))
    with pytest.raises(ValueError, match="rotational PTO coordinate must belong"):
        wec.rotational_pto("foreign", other.coordinate("pitch", foreign.move("pitch")))

    assert wec.bodies == [body]
    assert wec.coordinates == [coordinate]
    assert wec.rotational_ptos == []
    wec.pto("damper", body.at(0, 0, 0), WorldPoint(0, 0, 10), damping=500)
    case = wec.to_case(NoWave(), dt=0.1, end_time=0.2)
    assert [item["name"] for item in case["ptos"]] == ["damper"]


@pytest.mark.parametrize("settings,error", [
    ({"control": DeclutchingControl(100, 0.2), "stiffness": 10},
     "sampled PTO control sets its own gain"),
    ({"direct_drive": SimpleDirectDrive(1, 1, 1, 1, 1, 1, 1, 1),
      "damping": 10}, "direct drive supplies its own PTO force"),
    ({"linear_generator": DirectLinearGenerator(1, -1, 1, 1, 1, -1),
      "damping": 10}, "linear generator supplies its own PTO force"),
    ({"velocity_control_gain": -1}, "velocity_control_gain must be finite"),
    ({"position_control_gain": float("nan")},
     "position_control_gain must be finite"),
])
def test_pto_rejects_conflicting_or_invalid_force_laws_atomically(settings, error):
    wec, body, _ = _heave_wec()
    endpoints = (body.at(0, 0, 0), WorldPoint(0, 0, 10))
    with pytest.raises(ValueError, match=error):
        wec.pto("main", *endpoints, **settings)

    assert wec.ptos == []
    wec.pto("main", *endpoints, damping=500)
    with pytest.raises(ValueError, match="PTO name already exists"):
        wec.pto("main", *endpoints, damping=500)
    assert len(wec.ptos) == 1


def test_floating_joint_rejects_misordered_bodies_and_off_axis_location():
    wec = WEC("RM3")
    float_body = wec.body("float", HYDRO, inertia=(0, 1, 0))
    spar_body = wec.body("spar", HYDRO, inertia=(0, 1, 0))

    with pytest.raises(ValueError, match="exactly two ordered bodies"):
        wec.floating_joint(spar_body, float_body)
    with pytest.raises(ValueError, match="world point on the z axis"):
        wec.floating_joint(float_body, spar_body,
                           location=WorldPoint(1, 0, 0))
    with pytest.raises(ValueError, match="added-mass scheme"):
        wec.floating_joint(float_body, spar_body, added_mass_scheme="unknown")
    assert wec._floating_joint is None

    wec.floating_joint(float_body, spar_body, pto_name="relative_heave")
    case = wec.to_case(NoWave(), dt=0.1, end_time=0.2)
    assert case["constraint"]["kind"] == "floating_joint"
    assert case["pto"]["kind"] == "relative_heave"

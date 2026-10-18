"""Invalid public configurations must leave the WEC builder usable."""

from pathlib import Path

import numpy as np
import pytest

from wecsim import (
    DeclutchingControl,
    DirectLinearGenerator,
    HydroState,
    NoWave,
    PMWave,
    RegularWave,
    SimpleDirectDrive,
    WEC,
    WorldPoint,
)


HYDRO = Path(__file__).resolve().parents[1] / "examples/data/rm3.h5"
OSWEC = Path(__file__).resolve().parents[1] / "examples/data/oswec.h5"
GBM = Path(__file__).resolve().parents[1] / "tests/fixtures/matlab/gbm_barge.h5"
SPHERE = Path(__file__).resolve().parents[1] / "examples/data/sphere.h5"


def _heave_wec():
    wec = WEC("configured heave")
    body = wec.body("float", HYDRO)
    coordinate = wec.coordinate("heave", body.move("heave"))
    return wec, body, coordinate


def test_invalid_case_names_fail_before_changing_public_builder():
    wec = WEC("name recovery")
    states = (HydroState(HYDRO), HydroState(HYDRO))
    with pytest.raises(ValueError, match="body name"):
        wec.body("", HYDRO)
    with pytest.raises(ValueError, match="body name"):
        wec.fixed_body("", center_gravity=(0, 0, 0))
    with pytest.raises(ValueError, match="body name"):
        wec.variable_body("", states, switch_times=(.1,))
    assert not wec.bodies

    body = wec.body("float", HYDRO)
    with pytest.raises(ValueError, match="coordinate name"):
        wec.coordinate("bad name", body.move("heave"))
    assert not wec.coordinates
    coordinate = wec.coordinate("heave", body.move("heave"))
    endpoints = (WorldPoint(0, 0, -5), body.at(0, 0, 0))
    for invalid in ("bad name", "coordinate1"):
        with pytest.raises(ValueError, match="PTO name"):
            wec.pto(invalid, *endpoints, damping=500)
        with pytest.raises(ValueError, match="PTO name"):
            wec.rotational_pto(invalid, coordinate)
    assert not wec.ptos and not wec.rotational_ptos

    wec.pto("main", *endpoints, damping=500)
    result = wec.run(NoWave(), dt=.1, end_time=.2,
                     initial_coordinate={"heave": .1})
    assert result.ptos["main"].force.shape == (3,)


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


def test_variable_body_rejects_incomplete_state_clock_then_runs():
    wec = WEC("Sphere state clock recovery")
    states = (HydroState(SPHERE), HydroState(SPHERE))
    for configured_states, switch_times in (
        (states[:1], ()),
        (states, ()),
        (states, (.1, .2)),
    ):
        with pytest.raises(ValueError, match="states and one fewer switch times"):
            wec.variable_body("sphere", configured_states,
                              switch_times=switch_times)
        assert not wec.bodies

    body = wec.variable_body("sphere", states, switch_times=(.1,))
    wec.coordinate("heave", body.move("heave"))
    wec.pto("damper", WorldPoint(0, 0, 2), body.at(0, 0, 0),
            axis=(0, 0, 1), damping=200_000)
    result = wec.run(RegularWave(0, 8), dt=.01, end_time=.2, ramp_time=0,
                     initial_coordinate={"heave": .5})
    np.testing.assert_array_equal(
        dict(result.raw.extra_outputs)["body1_hydro_state"][10:], 2,
    )
    assert np.ptp(result.bodies["sphere"].position[:, 2]) > .01


def test_morison_attachment_recovery_prevents_silent_zero_load():
    wec = WEC("fixed Morison recovery")
    monopile = wec.fixed_body("monopile", center_gravity=(0, 0, -15))
    other = WEC("other")
    foreign = other.fixed_body("foreign", center_gravity=(0, 0, -15))
    parameters = dict(
        drag_coefficient=(1, 1, 1), added_mass_coefficient=(1, 1, 1),
        area=(100, 100, 100), volume=1000,
    )
    for body, point, axis, message in (
        (foreign, foreign.at(0, 0, 10), None, "must belong"),
        (monopile, foreign.at(0, 0, 10), None, "must belong"),
        (monopile, monopile.at(0, 0, 10), (0, 0, 0), "finite nonzero"),
    ):
        with pytest.raises(ValueError, match=message):
            wec.morison_element(body, point=point, axis=axis, **parameters)
        assert not wec._morison_elements

    wave = PMWave(2, 5, frequency_count=8, frequency_range=(.1, 2),
                  water_depth=30)
    with pytest.raises(ValueError, match="at least one element"):
        wec.run(wave, dt=.1, end_time=.2, ramp_time=0)
    wec.morison_element(monopile, point=monopile.at(0, 0, 10),
                        **parameters)
    result = wec.run(wave, dt=.1, end_time=.2, ramp_time=0)
    assert len(result.time) == 3
    assert np.max(np.abs(result.body_forces["monopile"])) > 0


def test_rotational_pto_rejects_translation_before_reserving_its_name():
    wec, body, heave = _heave_wec()
    with pytest.raises(ValueError, match="requires rotational motion"):
        wec.rotational_pto("hinge", heave, damping=500)

    assert wec.rotational_ptos == []
    pitch = wec.coordinate("pitch", body.move("pitch"))
    wec.rotational_pto("hinge", pitch, damping=500)
    case = wec.to_case(NoWave(), dt=.1, end_time=.2)
    assert case["ptos"] == [{
        "name": "hinge", "kind": "coordinate_torque", "coordinate": "pitch",
        "damping": 500, "stiffness": 0.0, "equilibrium_position": 0.0,
    }]


@pytest.mark.parametrize("damping", [True, np.bool_(True), -1, np.nan])
def test_rotational_pto_rejects_invalid_damping_before_reserving_name(damping):
    wec, body, _ = _heave_wec()
    pitch = wec.coordinate("pitch", body.move("pitch"))
    with pytest.raises(ValueError, match="finite nonnegative damping"):
        wec.rotational_pto("hinge", pitch, damping=damping)
    assert wec.rotational_ptos == []
    wec.rotational_pto("hinge", pitch, damping=500)
    assert wec.rotational_ptos[0].damping == 500


def test_coordinate_rejects_invalid_motion_before_reserving_its_name():
    wec, body, heave = _heave_wec()
    invalid = (
        ((body.move("diagonal"),), "unsupported DOF"),
        ((body.move("pitch", scale=0),), "finite nonzero scale"),
        ((body.move("pitch", scale=np.nan),), "finite nonzero scale"),
        ((body.move("pitch", scale=True),), "finite nonzero scale"),
        ((body.move("pitch"), body.move("pitch")), "repeats one body motion"),
    )
    for motions, message in invalid:
        with pytest.raises(ValueError, match=message):
            wec.coordinate("pitch", *motions)
        assert wec.coordinates == [heave]

    wec.coordinate("pitch", body.move("pitch"))
    case = wec.to_case(NoWave(), dt=.1, end_time=.2)
    assert [item["name"] for item in case["constraint"]["coordinates"]] == [
        "heave", "pitch",
    ]


@pytest.mark.parametrize("settings,error", [
    ({"control": DeclutchingControl(100, 0.2), "stiffness": 10},
     "sampled PTO control sets its own gain"),
    ({"direct_drive": SimpleDirectDrive(1, 1, 1, 1, 1, 1, 1, 1),
      "damping": 10}, "direct drive supplies its own PTO force"),
    ({"linear_generator": DirectLinearGenerator(1, -1, 1, 1, 1, -1),
      "damping": 10}, "linear generator supplies its own PTO force"),
    ({"velocity_control_gain": -1}, "velocity_control_gain must be finite"),
    ({"velocity_control_gain": np.bool_(True)},
     "velocity_control_gain must be finite"),
    ({"position_control_gain": float("nan")},
     "position_control_gain must be finite"),
    ({"position_control_gain": np.bool_(True)},
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


def test_linear_pto_rejects_invalid_attachment_and_force_before_reserving_name():
    wec, body, _ = _heave_wec()
    from_point = body.at(0, 0, 0)
    to_point = WorldPoint(0, 0, 10)
    invalid = (
        ((body.at(np.nan, 0, 0), to_point), {}, "finite three-vector"),
        ((from_point, WorldPoint(0, 0, np.inf)), {}, "finite three-vector"),
        ((from_point, to_point), {"axis": (0, 0, 0)}, "finite nonzero axis"),
        ((from_point, to_point), {"axis": (np.nan, 0, 1)}, "finite nonzero axis"),
        ((from_point, to_point), {"damping": -1}, "nonnegative damping"),
        ((from_point, to_point), {"stiffness": np.inf}, "stiffness must be finite"),
        ((from_point, to_point), {"pretension": 1}, "needs positive stiffness"),
    )
    for endpoints, settings, message in invalid:
        with pytest.raises(ValueError, match=message):
            wec.pto("main", *endpoints, **{"damping": 500, **settings})
        assert wec.ptos == []

    wec.pto("main", from_point, to_point, damping=500)
    result = wec.run(NoWave(), dt=.1, end_time=.2,
                     initial_coordinate={"heave": .1})
    assert result.ptos["main"].force.shape == (3,)


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


def test_floating_joint_rejects_invalid_solver_and_pto_without_committing_layout():
    wec = WEC("RM3 layout recovery")
    float_body = wec.body("float", HYDRO, inertia=(0, 21_306_090.66, 0))
    spar = wec.body("spar", HYDRO, inertia=(0, 94_407_091.24, 0))
    bad_settings = (
        ({"pto_name": ""}, ValueError, "PTO name"),
        ({"hard_stops": object()}, TypeError, "LinearHardStops"),
        ({"radiation_method": "unknown"}, ValueError, "radiation method"),
        ({"moordyn_point": spar.at(0, 0, 21.5)}, ValueError,
         "needs a MoorDyn session"),
    )
    for settings, error_type, message in bad_settings:
        with pytest.raises(error_type, match=message):
            wec.floating_joint(float_body, spar, **settings)
        assert wec._floating_joint is None

    wec.floating_joint(float_body, spar, damping=1_200_000)
    response = wec.run(RegularWave(0, 8), dt=.1, end_time=.2)
    assert response.case["constraint"]["kind"] == "floating_joint"
    assert response.ptos["relative_heave"].force.shape == (3,)
    np.testing.assert_allclose(
        response.ptos["relative_heave"].force,
        -1_200_000 * response.ptos["relative_heave"].velocity,
        rtol=0, atol=1e-9,
    )


def test_fixed_hinge_rejects_invalid_attachment_and_recovers():
    wec = WEC("OSWEC layout recovery")
    flap = wec.body("flap", OSWEC, mass=127_000,
                    inertia=(0, 1.85e6, 0))
    base = wec.fixed_body("base", center_gravity=(0, 0, -10.9), mass=999)
    bad_settings = (
        ({"location": WorldPoint(1, 0, -10)}, "z axis"),
        ({"pto_location": WorldPoint(0, 1, -8.9)}, "z axis"),
        ({"pto_name": ""}, "PTO name"),
        ({"damping": -1}, "nonnegative damping"),
        ({"added_mass_scheme": "unknown"}, "added-mass scheme"),
    )
    for settings, message in bad_settings:
        with pytest.raises(ValueError, match=message):
            wec.fixed_hinge(flap, base, **settings)
        assert wec._fixed_hinge is None

    wec.fixed_hinge(flap, base, stiffness=500, pretension=-50)
    response = wec.run(RegularWave(.1, 8), dt=.1, end_time=.2)
    assert response.case["pto"]["equilibrium_position"] == .1
    np.testing.assert_allclose(response.ptos["hinge"].force,
                               -500 * (response.ptos["hinge"].stroke - .1),
                               rtol=0, atol=1e-10)


def test_rejected_passive_yaw_bank_does_not_reserve_body_name():
    wec = WEC("heading bank recovery")
    for settings in (
        {"passive_yaw_threshold": 1},
        {"passive_yaw": True, "yaw_heading_bank": [0, 10, 10]},
        {"passive_yaw": True, "yaw_heading_bank": [-20, np.nan, 20]},
    ):
        with pytest.raises(ValueError, match="passive_yaw_threshold|yaw_heading_bank"):
            wec.body("flap", OSWEC, **settings)
        assert wec.bodies == []

    flap = wec.body("flap", OSWEC, passive_yaw=True,
                    yaw_heading_bank=[-20, 0, 20])
    assert wec.bodies == [flap]
    assert flap.yaw_heading_bank == (-20.0, 0.0, 20.0)


def test_rejected_gbm_force_settings_do_not_select_layout():
    wec = WEC("GBM recovery")
    barge = wec.body("barge", GBM, inertia=(6.667e7, 2.167e9, 2.167e9))
    for settings, error_type, message in (
        ({"orifice": object()}, TypeError, "OrificePTO"),
        ({"orifice_force_path": "published_owc"}, ValueError, "needs an orifice"),
        ({"heave_drag_area": 1}, ValueError, "need an orifice"),
        ({"heave_linear_damping": np.nan}, ValueError, "nonnegative and finite"),
        ({"pto_name": ""}, ValueError, "nonempty"),
    ):
        with pytest.raises(error_type, match=message):
            wec.floating_gbm(barge, **settings)

    wec.floating_gbm(barge)
    response = wec.run(RegularWave(0, 8), dt=.1, end_time=.2, ramp_time=0)
    assert response.case["constraint"]["kind"] == "floating_gbm"
    assert np.isfinite(response.bodies["barge"].position).all()


def test_rejected_physical_body_arrays_do_not_reserve_the_body_name():
    wec = WEC("body input recovery")
    valid_drag = (1, 1, 1, 0, 0, 0)
    invalid = (
        ({"hydrostatic_stiffness": "not a matrix"}, "finite 6x6 matrix"),
        ({"linear_damping": "not a matrix"}, "finite 6x6 matrix"),
        ({"quadratic_drag": "not a matrix"}, "finite 6x6 matrix"),
        ({"quadratic_drag_coefficient": ("bad",) * 6,
          "quadratic_drag_area": valid_drag}, "six finite nonnegative values"),
        ({"quadratic_drag_coefficient": valid_drag,
          "quadratic_drag_area": ("bad",) * 6}, "six finite nonnegative values"),
    )
    for settings, message in invalid:
        with pytest.raises(ValueError, match=message):
            wec.body("sphere", SPHERE, **settings)
        assert not wec.bodies

    damping = np.zeros((6, 6))
    damping[2, 2] = 100_000
    body = wec.body(
        "sphere", SPHERE, linear_damping=damping,
        quadratic_drag_coefficient=valid_drag,
        quadratic_drag_area=(100, 100, 100, 0, 0, 0),
    )
    wec.coordinate("heave", body.move("heave"))
    result = wec.run(NoWave(), dt=.05, end_time=.2,
                     radiation_memory=.1, initial_coordinate={"heave": .1})
    assert wec.bodies == [body]
    assert np.isfinite(result.bodies["sphere"].position).all()

    reference = WEC("undamped reference")
    plain = reference.body("sphere", SPHERE)
    reference.coordinate("heave", plain.move("heave"))
    undamped = reference.run(NoWave(), dt=.05, end_time=.2,
                             radiation_memory=.1,
                             initial_coordinate={"heave": .1})
    # The accepted body settings must affect dynamics, not merely serialize.
    assert (abs(result.bodies["sphere"].velocity[-1, 2])
            < abs(undamped.bodies["sphere"].velocity[-1, 2]))


def test_rejected_pto_force_settings_do_not_reserve_its_name():
    wec, body, _ = _heave_wec()
    endpoints = (WorldPoint(0, 0, -5), body.at(0, 0, 0))
    invalid = (
        ({"axis": ("bad", 0, 1)}, ValueError, "finite nonzero axis"),
        ({"axis": (0, 0, 0)}, ValueError, "finite nonzero axis"),
        ({"velocity_control_gain": "bad"}, ValueError,
         "velocity_control_gain must be finite"),
        ({"position_control_gain": "bad"}, ValueError,
         "position_control_gain must be finite"),
        ({"control": object()}, TypeError, "supported control law"),
        ({"direct_drive": object()}, TypeError, "SimpleDirectDrive"),
        ({"linear_generator": object()}, TypeError, "DirectLinearGenerator"),
    )
    for settings, error_type, message in invalid:
        with pytest.raises(error_type, match=message):
            wec.pto("main", *endpoints, **settings)
        assert not wec.ptos

    wec.pto("main", *endpoints, damping=500)
    result = wec.run(NoWave(), dt=.05, end_time=.2,
                     radiation_memory=.1, initial_coordinate={"heave": .1})
    assert np.isfinite(result.ptos["main"].force).all()

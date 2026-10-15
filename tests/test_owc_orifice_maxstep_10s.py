"""Separate OWC source solver-step sensitivity from piston-force routing."""

from hashlib import sha256
from pathlib import Path

import numpy as np
from scipy.spatial.transform import Rotation

from wecsim import OrificePTO, PMWave, WEC


FIXTURES = Path(__file__).parent / "fixtures/matlab"
HYDRO = FIXTURES / "owc_orifice.h5"
SOURCE = FIXTURES / "owc_orifice_maxstep_10s/source.npz"
PARAMETERS = FIXTURES / "owc_orifice_prefix.npz"


def _maximum_difference(actual, expected):
    assert actual.shape == expected.shape
    assert np.isfinite(actual).all() and np.isfinite(expected).all()
    return float(np.max(np.abs(actual - expected)))


def _public_result(force_path):
    with np.load(PARAMETERS, allow_pickle=False) as saved:
        orifice = OrificePTO(*saved["parameters"])
    wec = WEC("OWC internal-step audit")
    body = wec.body("OWC", HYDRO, inertia=(99.28, 11.04, 99.2))
    wec.floating_gbm(
        body, orifice=orifice, orifice_force_path=force_path,
        heave_linear_damping=100, mode_linear_damping=100,
        heave_drag_cd=1.2, heave_drag_area=8,
        pitch_drag_cd=1.2, pitch_drag_area=8,
    )
    return wec.run(
        PMWave(height=1, period=4, seed=128, phase_generator="matlab"),
        dt=.005, end_time=10, ramp_time=10, radiation_memory=15,
    )


def test_halving_only_matlab_maxstep_changes_owc_motion():
    assert sha256(SOURCE.read_bytes()).hexdigest() == (
        "d140edad47e60f230a72d7ba3cfe3628aa604bdddcdb99671b4aeb5587d03d35"
    )
    with np.load(SOURCE, allow_pickle=False) as saved:
        nominal = saved["body_nominal"]
        half = saved["body_half"]
        flex_nominal = saved["flex_nominal"]
        flex_half = saved["flex_half"]
        assert nominal.shape == half.shape == (2001, 13)
        assert flex_nominal.shape == flex_half.shape == (2001, 10)
        np.testing.assert_array_equal(saved["components_nominal"],
                                      saved["components_half"])
        np.testing.assert_array_equal(saved["wave_nominal"],
                                      saved["wave_half"])
        np.testing.assert_allclose(nominal[:, 0], np.arange(2001) * .005,
                                   rtol=0, atol=1e-12)
        np.testing.assert_array_equal(nominal[:, 0], half[:, 0])

    # The 0.1 mm and 0.1 mm/s source-stability limits were declared before
    # the local run. Euler roll/yaw can switch branches near pitch=pi/2, so
    # compare the physical relative orientation instead of raw Euler angles.
    nominal_rotation = Rotation.from_euler("xyz", nominal[:, 4:7])
    half_rotation = Rotation.from_euler("xyz", half[:, 4:7])
    angle_change = float(np.max((nominal_rotation.inv() * half_rotation).magnitude()))
    assert _maximum_difference(nominal[:, 1], half[:, 1]) > 1e-4
    assert _maximum_difference(nominal[:, 3], half[:, 3]) > 1e-4
    assert angle_change > 1e-4
    assert _maximum_difference(nominal[:, 7], half[:, 7]) > 1e-4
    assert _maximum_difference(flex_nominal[:, 0], flex_half[:, 0]) > 1e-4
    assert _maximum_difference(flex_nominal[:, 1], flex_half[:, 1]) > 1e-4


def test_public_owc_force_routes_against_halved_source_prefix():
    assert sha256(HYDRO.read_bytes()).hexdigest() == (
        "fe702ecf99f10d4cd0bba5d858b588dd37709cd45a5c3727356849e8dc40d676"
    )
    with np.load(SOURCE, allow_pickle=False) as saved:
        source_body = saved["body_half"]
        source_flex = saved["flex_half"]
        source_wave = saved["wave_half"]
    source_pitch = np.r_[0, np.cumsum(
        (source_body[1:, 11] + source_body[:-1, 11]) * .005 / 2,
    )]

    source_route = _public_result("published_owc")
    assert source_route.case["constraint"]["orifice_force_path"] == "published_owc"
    assert _maximum_difference(source_route.wave_elevation, source_wave[:, 1]) < 1e-11
    source_body_result = source_route.bodies["OWC"]
    assert _maximum_difference(source_body_result.position[:, 0], source_body[:, 1]) < .008
    assert _maximum_difference(source_body_result.position[:, 2], source_body[:, 3]) < .0015
    assert _maximum_difference(source_route.coordinates["pitch_unwrapped"].position,
                               source_pitch) < .007
    assert _maximum_difference(source_route.flexible_modes["OWC"].position[:, 0],
                               source_flex[:, 0]) < .0003

    # This is a diagnostic distinction, not a claim that the source route is
    # physically preferable or that the published 130 s case is paired.
    coupled = _public_result("coupled")
    assert coupled.case["constraint"]["orifice_force_path"] == "coupled"
    assert _maximum_difference(coupled.bodies["OWC"].position[:, 2],
                               source_body[:, 3]) > .05
    assert _maximum_difference(coupled.flexible_modes["OWC"].position[:, 0],
                               source_flex[:, 0]) > .015

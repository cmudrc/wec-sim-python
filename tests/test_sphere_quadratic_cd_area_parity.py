"""Gate WEC-Sim's density-scaled body-drag coefficient and area input."""

from pathlib import Path

import numpy as np
import pytest

from wecsim import NoWave, WEC
from wecsim.dynamics.caseDynamics import run_case


ROOT = Path(__file__).resolve().parents[1]
HYDRO = ROOT / "examples/data/sphere.h5"
SOURCE = ROOT / "tests/fixtures/matlab/sphere_quadratic_drag/source_cd_area.npz"
MATRIX_SOURCE = ROOT / "tests/fixtures/matlab/sphere_quadratic_drag/source.npz"
COEFFICIENT = (0, 0, 1, 0, 0, 0)
AREA = (0, 0, 500, 0, 0, 0)


def _heaving_sphere(**drag):
    wec = WEC("Sphere coefficient-area body drag")
    sphere = wec.body("sphere", HYDRO, **drag)
    wec.coordinate("heave", sphere.move("heave"))
    return wec


def test_coefficient_area_matches_pinned_matlab_sphere_free_decay():
    with np.load(SOURCE, allow_pickle=False) as source:
        time = source["time"].copy()
        position = source["position"].copy()
        velocity = source["velocity"].copy()
        force = source["viscous_force"].copy()
        coefficient = source["coefficient"].copy()
        area = source["area"].copy()
        rho = float(source["rho"])
    with np.load(MATRIX_SOURCE, allow_pickle=False) as matrix:
        np.testing.assert_allclose(time, matrix["time"], rtol=0, atol=1e-12)
        # The two distinct MATLAB input forms produce the same 40 s motion.
        np.testing.assert_allclose(position, matrix["position"], rtol=0, atol=1e-9)
        np.testing.assert_allclose(velocity, matrix["velocity"], rtol=0, atol=1e-9)
    np.testing.assert_array_equal(coefficient, COEFFICIENT)
    np.testing.assert_array_equal(area, AREA)
    assert rho == 1000
    diagonal = 0.5 * rho * coefficient * area
    assert np.max(np.abs(force - velocity * np.abs(velocity) * diagonal)) < 1e-5

    wec = _heaving_sphere(quadratic_drag_coefficient=COEFFICIENT,
                         quadratic_drag_area=AREA)
    result = wec.run(NoWave(), dt=.01, end_time=40, rho=rho,
                     radiation_memory=15, initial_coordinate={"heave": 1})
    np.testing.assert_allclose(result.time, time, rtol=0, atol=1e-12)
    assert result.case["bodies"][0]["quadratic_drag_coefficient"] == list(COEFFICIENT)
    assert result.case["bodies"][0]["quadratic_drag_area"] == list(AREA)
    actual = result.bodies["sphere"]
    assert np.max(np.abs(actual.position[:, 2] - position[:, 2])) < 1e-4
    assert np.max(np.abs(actual.velocity[:, 2] - velocity[:, 2])) < 1e-4
    for dof in (0, 1, 3, 4, 5):
        assert np.max(np.abs(actual.position[:, dof] - position[:, dof])) < 1e-10
        assert np.max(np.abs(actual.velocity[:, dof] - velocity[:, dof])) < 1e-10


def test_coefficient_area_scales_with_run_density():
    rho = 1025
    coefficient = _heaving_sphere(quadratic_drag_coefficient=COEFFICIENT,
                                 quadratic_drag_area=AREA)
    matrix = np.diag(0.5 * rho * np.asarray(COEFFICIENT) * np.asarray(AREA))
    explicit = _heaving_sphere(quadratic_drag=matrix)
    settings = dict(dt=.01, end_time=2, rho=rho, radiation_memory=2,
                    initial_coordinate={"heave": 1})
    first = coefficient.run(NoWave(), **settings).bodies["sphere"]
    second = explicit.run(NoWave(), **settings).bodies["sphere"]
    np.testing.assert_allclose(first.position, second.position, rtol=0, atol=1e-11)
    np.testing.assert_allclose(first.velocity, second.velocity, rtol=0, atol=1e-11)


@pytest.mark.parametrize("drag, message", [
    ({"quadratic_drag_coefficient": COEFFICIENT}, "supplied together"),
    ({"quadratic_drag_area": AREA}, "supplied together"),
    ({"quadratic_drag_coefficient": (0, 0, True, 0, 0, 0),
      "quadratic_drag_area": AREA}, "six finite nonnegative"),
    ({"quadratic_drag_coefficient": COEFFICIENT,
      "quadratic_drag_area": (0, 0, -1, 0, 0, 0)}, "six finite nonnegative"),
    ({"quadratic_drag_coefficient": COEFFICIENT,
      "quadratic_drag_area": AREA, "quadratic_drag": np.zeros((6, 6))}, "either"),
])
def test_public_body_rejects_incomplete_or_ambiguous_drag(drag, message):
    wec = WEC("body drag validation")
    with pytest.raises(ValueError, match=message):
        wec.body("sphere", HYDRO, **drag)
    assert wec.bodies == []


def test_direct_case_rejects_incomplete_coefficient_area_pair():
    case = _heaving_sphere(quadratic_drag_coefficient=COEFFICIENT,
                          quadratic_drag_area=AREA).to_case(
                              NoWave(), dt=.01, end_time=.1, radiation_memory=.1,
                          )
    del case["bodies"][0]["quadratic_drag_area"]
    with pytest.raises(ValueError, match="supplied together"):
        run_case(case)


def test_direct_case_rejects_null_drag_vectors():
    case = _heaving_sphere(quadratic_drag_coefficient=COEFFICIENT,
                          quadratic_drag_area=AREA).to_case(
                              NoWave(), dt=.01, end_time=.1, radiation_memory=.1,
                          )
    case["bodies"][0]["quadratic_drag_coefficient"] = None
    case["bodies"][0]["quadratic_drag_area"] = None
    with pytest.raises(ValueError, match="six finite nonnegative"):
        run_case(case)


def test_specialized_joint_rejects_unmodeled_body_drag_vectors():
    wec = WEC("RM3 body drag scope")
    hydro = ROOT / "examples/data/rm3.h5"
    float_body = wec.body(
        "float", hydro, inertia=(0, 21_306_090.66, 0),
        quadratic_drag_coefficient=COEFFICIENT, quadratic_drag_area=AREA,
    )
    spar = wec.body("spar", hydro, inertia=(0, 94_407_091.24, 0))
    wec.floating_joint(float_body, spar)
    with pytest.raises(ValueError, match="without extra force models"):
        wec.to_case(NoWave(), dt=.1, end_time=.2)

"""Keep the published MATLAB Sphere free-decay cases paired in Python CI."""

from pathlib import Path

import numpy as np
import pytest

from wecsim import NoWave, WEC
from wecsim.dynamics.caseDynamics import run_case


ROOT = Path(__file__).resolve().parents[1]
HYDRO = ROOT / "examples/data/sphere.h5"
REFERENCE = ROOT / "tests/fixtures/matlab/sphere_free_decay.npz"
CASES = (("0m", 0), ("1m", 1), ("1m-ME", 1), ("3m", 3), ("5m", 5))


@pytest.fixture(scope="module")
def source():
    with np.load(REFERENCE, allow_pickle=False) as saved:
        arrays = {name: saved[name] for name in saved.files}
    np.testing.assert_array_equal(arrays["case_names"], [case for case, _ in CASES])
    assert arrays["time"].shape == (4001,)
    assert arrays["body_position"].shape == (5, 4001, 6)
    assert arrays["body_velocity"].shape == (5, 4001, 6)
    assert arrays["body_total_heave_force"].shape == (5, 4001)
    return arrays


def _max_error(actual, expected, limit, label):
    assert actual.shape == expected.shape, f"{label}: mismatched shape"
    assert np.isfinite(actual).all() and np.isfinite(expected).all(), label
    error = float(np.max(np.abs(actual - expected)))
    assert error < limit, f"{label}: {error:.6g} exceeds {limit:.6g}"


@pytest.mark.parametrize("index,case,displacement", [
    (index, case, displacement)
    for index, (case, displacement) in enumerate(CASES)
])
def test_published_sphere_free_decay(source, index, case, displacement):
    response = run_case({
        "simulation": {"dt": 0.01, "end_time": 40, "radiation_memory": 15},
        "wave": {"type": "none"},
        "bodies": [{"hydro_file": str(HYDRO), "mass": "equilibrium"}],
        "constraint": {"kind": "heave", "initial_displacement": displacement},
    })
    np.testing.assert_allclose(response.time, source["time"], rtol=0,
                               atol=1e-10)
    scale = max(1, displacement)
    for dof in range(6):
        position_limit = 1e-10 if dof != 2 or displacement == 0 else 1.2e-4 * scale
        velocity_limit = 1e-10 if dof != 2 or displacement == 0 else 1.6e-4 * scale
        _max_error(response.body_position[:, 0, dof],
                   source["body_position"][index, :, dof],
                   position_limit, f"Sphere {case} position DOF {dof}")
        _max_error(response.body_velocity[:, 0, dof],
                   source["body_velocity"][index, :, dof],
                   velocity_limit, f"Sphere {case} speed DOF {dof}")
    force_limit = 1e-10 if displacement == 0 else 100 * scale
    _max_error(response.total_heave_force,
               source["body_total_heave_force"][index],
               force_limit, f"Sphere {case} total heave force")


def test_public_wec_configures_published_five_metre_free_decay(source):
    wec = WEC("Sphere free decay")
    sphere = wec.body("sphere", HYDRO)
    wec.coordinate("heave", sphere.move("heave"))
    result = wec.run(NoWave(), dt=0.01, end_time=40, radiation_memory=15,
                     initial_coordinate={"heave": 5})
    np.testing.assert_allclose(result.time, source["time"], rtol=0,
                               atol=1e-10)
    index = len(CASES) - 1
    for dof in range(6):
        limit = 0.001 if dof == 2 else 1e-10
        _max_error(result.bodies["sphere"].position[:, dof],
                   source["body_position"][index, :, dof], limit,
                   f"public Sphere position DOF {dof}")
        _max_error(result.bodies["sphere"].velocity[:, dof],
                   source["body_velocity"][index, :, dof], limit,
                   f"public Sphere speed DOF {dof}")

"""Exercise the case input boundary independently of MATLAB availability."""

import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest
from scipy.io import savemat

from wecsim.dynamics.caseDynamics import run_case


ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / "examples/rm3.json"


def test_example_case_runs_and_records_provenance(tmp_path):
    case = json.loads(EXAMPLE.read_text(encoding="utf-8"))
    case["simulation"]["end_time"] = 2
    case_file = tmp_path / "case.json"
    for body in case["bodies"]:
        body["hydro_file"] = str((EXAMPLE.parent / body["hydro_file"]).resolve())
    case_file.write_text(json.dumps(case), encoding="utf-8")
    output = tmp_path / "motion.csv"
    subprocess.run(
        [sys.executable, "-m", "wecsim",
         str(case_file), "--output", str(output)],
        cwd=ROOT, check=True, capture_output=True, text=True,
    )
    values = np.loadtxt(output, delimiter=",", skiprows=1)
    metadata = json.loads(output.with_suffix(".json").read_text(encoding="utf-8"))
    assert values.shape == (21, 27)
    assert np.isfinite(values).all()
    assert metadata["time_steps"] == 21
    assert len(metadata["hydro_files"]) == 1
    assert len(metadata["hydro_files"][0]["sha256"]) == 64


def test_cli_keeps_case_when_csv_has_same_stem(tmp_path):
    case = json.loads(EXAMPLE.read_text(encoding="utf-8"))
    case["simulation"].update(end_time=0.2, ramp_time=0)
    for body in case["bodies"]:
        body["hydro_file"] = str((EXAMPLE.parent / body["hydro_file"]).resolve())
    case_file = tmp_path / "case.json"
    original = json.dumps(case)
    case_file.write_text(original, encoding="utf-8")
    output = tmp_path / "case.csv"

    subprocess.run(
        [sys.executable, "-m", "wecsim", str(case_file),
         "--output", str(output)],
        cwd=ROOT, check=True, capture_output=True, text=True,
    )

    assert case_file.read_text(encoding="utf-8") == original
    values = np.loadtxt(output, delimiter=",", skiprows=1)
    assert values.shape == (3, 27) and np.isfinite(values).all()
    metadata = json.loads(
        (tmp_path / "case.provenance.json").read_text(encoding="utf-8"),
    )
    assert metadata["case"] == case
    assert metadata["time_steps"] == 3
    rejected = subprocess.run(
        [sys.executable, "-m", "wecsim", str(case_file),
         "--output", str(case_file)],
        cwd=ROOT, capture_output=True, text=True,
    )
    assert rejected.returncode == 2
    assert "cannot overwrite the case input" in rejected.stderr
    assert case_file.read_text(encoding="utf-8") == original


def test_independent_multibody_radiation_maps_each_body_to_its_own_velocity():
    hydro = str((EXAMPLE.parent / json.loads(EXAMPLE.read_text())[
        "bodies"][0]["hydro_file"]).resolve())
    first = [[0, 0], [0, 0], [1, 0], [0, 0], [0, 0], [0, 0]]
    second = [[0, 0], [0, 0], [0, 1], [0, 0], [0, 0], [0, 0]]
    case = {
        "simulation": {"dt": 0.1, "end_time": 0.3,
                       "radiation_memory": 0.2},
        "wave": {"type": "none"},
        "bodies": [
            {"hydro_file": hydro, "coordinate_map": first},
            {"hydro_file": hydro, "coordinate_map": second},
        ],
        "constraint": {"kind": "linear_subspace",
                       "initial_coordinate": [1, 0]},
    }
    independent = run_case(case)
    one_body = {
        **case,
        "bodies": [{"hydro_file": hydro,
                    "coordinate_map": [[row[0]] for row in first]}],
        "constraint": {"kind": "linear_subspace",
                       "initial_coordinate": [1]},
    }
    reference = run_case(one_body)
    np.testing.assert_allclose(
        independent.body_position[:, 0, 2],
        reference.body_position[:, 0, 2], rtol=0, atol=1e-10,
    )
    np.testing.assert_allclose(
        independent.body_position[:, 1, 2],
        independent.body_position[0, 1, 2], rtol=0, atol=1e-10,
    )
    coupled = run_case({**case, "body_to_body": True})
    assert np.isfinite(coupled.body_position).all()


def test_scalar_pto_equilibrium_and_pretension_change_force():
    case = json.loads(EXAMPLE.read_text(encoding="utf-8"))
    case["simulation"]["end_time"] = 0.2
    case["pto"].update(damping=0, stiffness=1000,
                       equilibrium_position=0.5)
    shifted = run_case(case, base_dir=EXAMPLE.parent)
    np.testing.assert_allclose(shifted.pto_force[0], 500, rtol=0, atol=1e-12)
    case["pto"]["equilibrium_position"] = 0
    unshifted = run_case(case, base_dir=EXAMPLE.parent)
    assert np.max(np.abs(shifted.body_position - unshifted.body_position)) > 1e-9
    del case["pto"]["equilibrium_position"]
    case["pto"]["pretension"] = 100
    pretensioned = run_case(case, base_dir=EXAMPLE.parent)
    np.testing.assert_allclose(pretensioned.pto_force[0], -100, rtol=0,
                               atol=1e-12)
    case["pto"]["equilibrium_position"] = 0.5
    with pytest.raises(ValueError, match="either PTO equilibrium_position or pretension"):
        run_case(case, base_dir=EXAMPLE.parent)


def test_mapped_pto_equilibrium_changes_generalized_force():
    hydro = json.loads(EXAMPLE.read_text())["bodies"][0]["hydro_file"]
    case = {
        "simulation": {"dt": 0.1, "end_time": 0, "ramp_time": 0},
        "wave": {"type": "regular", "height": 0, "period": 8},
        "bodies": [{"hydro_file": hydro,
                    "coordinate_map": [[0], [0], [1], [0], [0], [0]]}],
        "constraint": {"kind": "linear_subspace"},
        "pto": {"kind": "linear", "stiffness_matrix": [[1000]],
                "equilibrium_coordinate": [0.5]},
    }
    response = run_case(case, base_dir=EXAMPLE.parent)
    np.testing.assert_allclose(response.pto_generalized_force, [[500]],
                               rtol=0, atol=1e-12)


@pytest.mark.parametrize("option", [
    "control_surface", "momentum_conservation", "pressure_integration",
])
def test_mean_drift_requires_the_selected_h5_coefficient(option):
    hydro = ROOT / "examples/data/rm3.h5"
    case = {
        "simulation": {"dt": 0.1, "end_time": 0,
                       "radiation_memory": 0.1},
        "wave": {"type": "regularCIC", "height": 0.1, "period": 8},
        "bodies": [{"hydro_file": str(hydro),
                    "mean_drift": option,
                    "coordinate_map": [[0], [0], [1], [0], [0], [0]]}],
        "constraint": {"kind": "linear_subspace"},
    }
    with pytest.raises(ValueError, match=f"lacks finite {option}"):
        run_case(case)


@pytest.mark.parametrize("change,explanation", [
    ({"constraint": {"kind": "free_six_dof"}}, "unsupported constraint layout"),
    ({"pto": {"kind": "pitch", "damping": 1}}, "relative_heave PTO"),
    ({"wave": {"type": "irregular"}}, "regular waves"),
])
def test_unsupported_physics_fails_explicitly(change, explanation):
    case = json.loads(EXAMPLE.read_text(encoding="utf-8"))
    case["simulation"]["end_time"] = 0
    case.update(change)
    with pytest.raises(ValueError, match=explanation):
        run_case(case, base_dir=EXAMPLE.parent)


def test_imported_wave_rejects_incomplete_time_coverage(tmp_path):
    case = json.loads(EXAMPLE.read_text(encoding="utf-8"))
    case["simulation"].update(dt=0.1, end_time=1, ramp_time=0)
    case["wave"] = {"type": "elevationImport", "file": "short.mat"}
    for body in case["bodies"]:
        body["hydro_file"] = str((EXAMPLE.parent / body["hydro_file"]).resolve())
    savemat(tmp_path / "short.mat", {"etaData": [[0, 0], [0.5, 0]]})
    with pytest.raises(ValueError, match="cover the full simulation time"):
        run_case(case, base_dir=tmp_path)

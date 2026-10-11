"""Check the public reference-case command and its reproducibility record."""

import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize("model,h5_file,columns,extra", [
    ("rm3", ROOT / "examples" / "data" / "rm3.h5", 14, []),
    ("rm3", ROOT / "examples" / "data" / "rm3.h5", 14, ["--b2b"]),
    ("oswec", ROOT / "tests" / "test_objects" / "test_bodyclass"
     / "testData" / "hydroData" / "oswec.h5", 7, []),
])
def test_reference_runner_exports_result_and_metadata(tmp_path, model, h5_file, columns, extra):
    output = tmp_path / f"{model}.csv"
    subprocess.run(
        [sys.executable, "-m", "wecsim.reference", model,
         "--h5", str(h5_file), "--output", str(output),
         "--end-time", "2", *extra],
        cwd=ROOT, check=True, capture_output=True, text=True,
    )
    values = np.loadtxt(output, delimiter=",", skiprows=1)
    metadata = json.loads(output.with_suffix(".json").read_text())
    assert values.shape == (21, columns)
    assert np.isfinite(values).all()
    assert metadata["model"] == model
    assert metadata["time_steps"] == 21
    assert len(metadata["h5_sha256"]) == 64
    assert len(metadata["csv_columns"]) == columns
    if model == "rm3":
        assert metadata["parameters"]["b2b"] == ("--b2b" in extra)


def test_reference_runner_keeps_hydrodynamic_input(tmp_path):
    h5_file = tmp_path / "rm3.h5"
    h5_file.write_bytes((ROOT / "examples/data/rm3.h5").read_bytes())
    original = h5_file.read_bytes()
    rejected = subprocess.run(
        [sys.executable, "-m", "wecsim.reference", "rm3",
         "--h5", str(h5_file), "--output", str(h5_file)],
        cwd=ROOT, capture_output=True, text=True,
    )
    assert rejected.returncode == 2
    assert "cannot overwrite the hydrodynamic input" in rejected.stderr
    assert h5_file.read_bytes() == original

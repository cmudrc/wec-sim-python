"""Exercise the installed command's CSV/provenance contract with a short case."""

import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from wecsim.interface.cli import main


ROOT = Path(__file__).resolve().parents[1]
HYDRO = ROOT / "examples/data/rm3.h5"
EXAMPLE = ROOT / "examples/configurable_rm3_pto.json"


def _case(tmp_path):
    case = json.loads(EXAMPLE.read_text())
    case["simulation"].update(end_time=.2, ramp_time=0)
    for body in case["bodies"]:
        body["hydro_file"] = str(HYDRO)
    path = tmp_path / "case.json"
    path.write_text(json.dumps(case))
    return path


def test_cli_exports_short_case_with_traceable_inputs(tmp_path, capsys):
    case = _case(tmp_path)
    output = tmp_path / "case.csv"
    main([str(case), "--output", str(output)])
    assert capsys.readouterr().out.strip() == str(output)
    header = output.read_text().splitlines()[0].split(",")
    data = np.loadtxt(output, delimiter=",", skiprows=1)
    assert data.shape == (3, len(header))
    assert header[:2] == ["time", "wave_elevation"]
    assert "pto_main_force" in header
    np.testing.assert_allclose(data[:, 0], [0, .1, .2])
    assert np.isfinite(data).all()

    # The default .json provenance name collides with case.json, so the
    # command must preserve its input and use the fallback name.
    assert json.loads(case.read_text())["simulation"]["end_time"] == .2
    metadata = json.loads((tmp_path / "case.provenance.json").read_text())
    assert metadata["case_sha256"] == hashlib.sha256(case.read_bytes()).hexdigest()
    assert metadata["hydro_files"] == [{
        "path": str(HYDRO),
        "sha256": hashlib.sha256(HYDRO.read_bytes()).hexdigest(),
    }]
    assert metadata["csv_columns"] == header
    assert metadata["time_steps"] == 3


def test_cli_refuses_to_replace_case_or_hydro_input(tmp_path):
    case = _case(tmp_path)
    with pytest.raises(SystemExit, match="2"):
        main([str(case), "--output", str(case)])
    with pytest.raises(SystemExit, match="2"):
        main([str(case), "--output", str(HYDRO)])

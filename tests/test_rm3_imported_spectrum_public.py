"""Full published RM3 imported-spectrum motion through the public WEC API."""

from hashlib import sha256
from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import trapezoid

from wecsim import ImportedSpectrumWave, WEC


ROOT = Path(__file__).resolve().parents[1]
HYDRO = ROOT / "examples/data/rm3.h5"
FIXTURE = ROOT / "tests/fixtures/matlab/rm3_mcr_seastate"
SOURCE = FIXTURE / "source.npz"


def _max_error(actual, expected, limit, label):
    assert actual.shape == expected.shape, label
    assert np.isfinite(actual).all() and np.isfinite(expected).all(), label
    difference = float(np.max(np.abs(actual - expected)))
    assert difference < limit, f"{label}: {difference:.6g} >= {limit:.6g}"


def test_archived_imported_spectrum_baseline_provenance():
    assert sha256(SOURCE.read_bytes()).hexdigest() == (
        "501145af44398328287f1e3d2fd15a8b0ee3c8e84c22f5f04685279dbaa8f6d7"
    )
    with np.load(SOURCE, allow_pickle=False) as saved:
        assert saved["summary"].shape == (3,)
        for index in (1, 2, 3):
            assert saved[f"case{index}"].shape == (4001, 18)


def _wec(index, *, source_timing):
    wec = WEC(f"RM3 imported spectrum {index}")
    float_body = wec.body("float", HYDRO, inertia=(0, 21_306_090.66, 0))
    spar = wec.body("spar", HYDRO, inertia=(0, 94_407_091.24, 0))
    settings = {"added_mass_scheme": "simulink_delay"} if source_timing else {}
    wec.floating_joint(float_body, spar, damping=1_200_000, **settings)
    return wec


@pytest.mark.parametrize("index", (1, 2, 3))
def test_public_floating_joint_reproduces_imported_spectrum_mcr(index):
    with np.load(SOURCE, allow_pickle=False) as saved:
        source = saved[f"case{index}"].copy()
        source_mean_power = float(saved["summary"][index - 1])
    spectrum = FIXTURE / f"spectrumData{index}.mat"
    result = _wec(index, source_timing=True).run(
        ImportedSpectrumWave(spectrum), dt=.1, end_time=400,
        ramp_time=100, radiation_memory=60,
    )
    assert result.case["wave"]["type"] == "spectrumImport"
    assert result.case["simulation"]["added_mass_scheme"] == "simulink_delay"
    assert result.raw.auxiliary_files == (spectrum.resolve(),)
    np.testing.assert_allclose(result.time, source[:, 0], rtol=0, atol=1e-10)
    _max_error(result.wave_elevation, source[:, 1], 1e-10,
               f"case {index} elevation")

    for name, offset in (("float", 2), ("spar", 8)):
        body = result.bodies[name]
        for dof, axis, position_limit, speed_limit in (
            (0, "surge", .020, .005),
            (2, "heave", .0015, .0013),
            (4, "pitch", .0009, .0004),
        ):
            column = {0: 0, 2: 1, 4: 2}[dof]
            _max_error(body.position[:, dof], source[:, offset + column],
                       position_limit, f"case {index} {name} {axis} position")
            _max_error(body.velocity[:, dof], source[:, offset + 3 + column],
                       speed_limit, f"case {index} {name} {axis} speed")

    pto = result.ptos["relative_heave"]
    for label, actual, expected, limit in (
        ("stroke", pto.stroke, source[:, 14], .0035),
        ("speed", pto.velocity, source[:, 15], .002),
        ("force", pto.force, source[:, 16], 2_500),
        ("power", pto.absorbed_power, -source[:, 17], 1_800),
    ):
        _max_error(actual, expected, limit, f"case {index} PTO {label}")
    source_work = -trapezoid(source[:, 17], x=source[:, 0])
    python_work = trapezoid(pto.absorbed_power, x=result.time)
    assert source_work > 0
    assert abs(python_work - source_work) / source_work < .001
    assert abs(np.mean(pto.absorbed_power[1999:]) + source_mean_power) < 75


@pytest.mark.parametrize("index", (1, 2, 3))
def test_physical_default_keeps_imported_spectrum_pto_work_dissipative(index):
    with np.load(SOURCE, allow_pickle=False) as saved:
        source = saved[f"case{index}"].copy()
    spectrum = FIXTURE / f"spectrumData{index}.mat"
    result = _wec(index, source_timing=False).run(
        ImportedSpectrumWave(spectrum), dt=.1, end_time=400,
        ramp_time=100, radiation_memory=60,
    )
    assert result.case["simulation"]["added_mass_scheme"] == "implicit"
    pto = result.ptos["relative_heave"]
    assert np.all(pto.absorbed_power >= 0)
    source_work = -trapezoid(source[:, 17], x=source[:, 0])
    physical_work = trapezoid(pto.absorbed_power, x=result.time)
    assert source_work > 0
    assert abs(physical_work - source_work) / source_work < .01

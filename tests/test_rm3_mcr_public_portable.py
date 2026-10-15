"""Portable public-WEC gates for all eight published RM3 MCR conditions."""

from hashlib import sha256
from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import trapezoid

from wecsim import RegularCICWave, WEC


ROOT = Path(__file__).resolve().parents[1]
HYDRO = ROOT / "examples/data/rm3.h5"
SOURCE = ROOT / "tests/fixtures/matlab/rm3_mcr_array/source.npz"


def _source(index):
    with np.load(SOURCE, allow_pickle=False) as saved:
        return saved["summary"][index - 1].copy(), saved[f"case{index}"].copy()


def _run(index, parameters, *, source_timing):
    height, period, damping = parameters[:3]
    wec = WEC(f"RM3 MCR condition {index}")
    float_body = wec.body("float", HYDRO, inertia=(0, 21_306_090.66, 0))
    spar = wec.body("spar", HYDRO, inertia=(0, 94_407_091.24, 0))
    kwargs = {"added_mass_scheme": "simulink_delay"} if source_timing else {}
    wec.floating_joint(float_body, spar, damping=float(damping),
                       radiation_method="convolution", **kwargs)
    return wec.run(RegularCICWave(float(height), float(period)),
                   dt=.1, end_time=400, ramp_time=100, radiation_memory=60)


def _max_error(actual, expected, limit, label):
    assert actual.shape == expected.shape, label
    assert np.isfinite(actual).all() and np.isfinite(expected).all(), label
    error = float(np.max(np.abs(actual - expected)))
    assert error < limit, f"{label}: {error:.6g} >= {limit:.6g}"


def test_archived_mcr_array_provenance_and_condition_order():
    assert sha256(SOURCE.read_bytes()).hexdigest() == (
        "2de50eb304b39c7fe37f08afc9d55168e9d461655ce3b29ae5bceb77f8a36ed1"
    )
    with np.load(SOURCE, allow_pickle=False) as saved:
        summary = saved["summary"]
        assert summary.shape == (8, 6)
        np.testing.assert_array_equal(summary[:, 0], [1.5, 1.5, 2.5, 2.5] * 2)
        np.testing.assert_array_equal(summary[:, 1], [6, 8, 6, 8] * 2)
        np.testing.assert_array_equal(summary[:, 2], [1_200_000] * 4 + [2_400_000] * 4)
        np.testing.assert_array_equal(summary[:, 3], 0)


@pytest.mark.parametrize("index", range(1, 9))
def test_public_source_timing_pairs_each_published_mcr_condition(index):
    parameters, source = _source(index)
    assert source.shape == (4001, 17)
    result = _run(index, parameters, source_timing=True)
    np.testing.assert_allclose(result.time, source[:, 0], rtol=0, atol=1e-10)
    assert result.case["simulation"]["added_mass_scheme"] == "simulink_delay"
    assert result.case["pto"]["damping"] == parameters[2]

    for name, offset in (("float", 1), ("spar", 7)):
        body = result.bodies[name]
        for dof, axis, position_limit, speed_limit in (
            (0, "surge", .09, .003),
            (2, "heave", .0006, .0005),
            (4, "pitch", .0003, .00012),
        ):
            column = {0: 0, 2: 1, 4: 2}[dof]
            _max_error(body.position[:, dof], source[:, offset + column],
                       position_limit, f"case {index} {name} {axis} position")
            _max_error(body.velocity[:, dof], source[:, offset + 3 + column],
                       speed_limit, f"case {index} {name} {axis} speed")

    pto = result.ptos["relative_heave"]
    for label, actual, expected, limit in (
        ("stroke", pto.stroke, source[:, 13], .0006),
        ("speed", pto.velocity, source[:, 14], .0004),
        ("force", pto.force, source[:, 15], 700),
        ("power", pto.absorbed_power, -source[:, 16], 700),
    ):
        _max_error(actual, expected, limit, f"case {index} PTO {label}")
    source_work = -trapezoid(source[:, 16], x=source[:, 0])
    python_work = trapezoid(pto.absorbed_power, x=result.time)
    assert source_work > 0
    assert abs(python_work - source_work) / source_work < .001
    assert abs(np.mean(pto.absorbed_power[1999:]) + parameters[4]) < 250


@pytest.mark.parametrize("index", range(1, 9))
def test_physical_default_keeps_mcr_pto_work_dissipative(index):
    parameters, source = _source(index)
    result = _run(index, parameters, source_timing=False)
    assert result.case["simulation"].get("added_mass_scheme", "implicit") == "implicit"
    pto = result.ptos["relative_heave"]
    assert np.all(pto.absorbed_power >= 0)
    source_work = -trapezoid(source[:, 16], x=source[:, 0])
    python_work = trapezoid(pto.absorbed_power, x=result.time)
    assert source_work > 0
    assert abs(python_work - source_work) / source_work < .01

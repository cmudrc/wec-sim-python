"""Separate RM3 PTO-extension ode4 integration from WEC-Sim sample clocks."""

import json
import os
from pathlib import Path

import numpy as np
import pytest

from wecsim import NoWave, WEC


ROOT = Path(__file__).resolve().parents[1]
RECORDS = os.environ.get("WEC_SIM_RM3_PTO_FIXEDSTEP_DIR")
SELECTED_CASE = os.environ.get("WEC_SIM_RM3_PTO_CASE")
HYDRO = ROOT / "examples/data/rm3_end_stops.h5"
ARCHIVE = ROOT / "tests/fixtures/matlab/rm3_pto_extension_step.npz"
FIXTURE = ROOT / "tests/fixtures/matlab/rm3_pto_extension_fixedstep_only/source.npz"


def _record(folder, name, mode):
    if folder is None:
        with np.load(FIXTURE, allow_pickle=False) as source:
            values = source[f"{name}_{mode}"]
            settings = source[f"{name}_{mode}_settings"].reshape(1, 5)
    else:
        values = np.loadtxt(folder / f"{mode}.csv", delimiter=",", ndmin=2)
        settings = np.loadtxt(folder / f"{mode}_settings.csv", delimiter=",", ndmin=2)
    assert values.shape == (301, 9) and np.isfinite(values).all()
    assert settings.shape == (1, 5) and np.isfinite(settings).all()
    np.testing.assert_allclose(values[:, 0], np.arange(301) * .1,
                               rtol=0, atol=1e-10)
    np.testing.assert_allclose(settings[0],
                               [.1, .1, .1, .1,
                                .1 if mode == "published" else .0125],
                               rtol=0, atol=1e-12)
    return values


def _max_gap(a, b):
    return float(np.max(np.abs(a - b)))


def _physical_run(name):
    wec = WEC(f"RM3 {name} PTO-extension fixed-step audit")
    float_body = wec.body("float", HYDRO, inertia=(0, 21_306_090.66, 0))
    spar = wec.body("spar", HYDRO, inertia=(0, 94_407_091.24, 0))
    wec.floating_joint(float_body, spar, added_mass_scheme="implicit")
    initial = {"float_heave": 5} if name == "float" else {"spar_heave": -5}
    return wec.run(NoWave(), dt=.0125, end_time=30,
                   radiation_memory=60, initial_coordinate=initial)


@pytest.mark.parametrize("name,column", (
    [(SELECTED_CASE, 1 if SELECTED_CASE == "float" else 3)]
    if SELECTED_CASE else [("float", 1), ("spar", 3)]
))
def test_fixedstep_alone_explains_most_published_extension_gap(name, column):
    folder = Path(RECORDS) / name if RECORDS else None
    published = _record(folder, name, "published")
    fixed = _record(folder, name, "fixedstep0125")
    with np.load(ARCHIVE, allow_pickle=False) as source:
        prior_published = source[f"{name}_0p1"]
        refined = source[f"{name}_0p0125"][::8]
    assert prior_published.shape == refined.shape == (301, 9)
    np.testing.assert_allclose(published, prior_published, rtol=0, atol=1e-7)
    np.testing.assert_allclose(refined[:, 0], published[:, 0], rtol=0, atol=1e-10)

    # The published PTO is unforced and the other body remains stationary.
    assert _max_gap(published[:, 7:9], 0) < 1e-7
    assert _max_gap(fixed[:, 7:9], 0) < 1e-7
    inactive = 3 if name == "float" else 1
    assert _max_gap(fixed[:, inactive], fixed[0, inactive]) < 1e-8
    assert _max_gap(fixed[:, inactive + 1], 0) < 1e-8

    metrics = {}
    for label, index in (("heave", column), ("speed", column + 1)):
        original_error = refined[:, index] - published[:, index]
        residual = refined[:, index] - fixed[:, index]
        squared_reduction = 1 - float(residual @ residual / (original_error @ original_error))
        metrics[f"source_published_to_refined_{label}"] = _max_gap(original_error, 0)
        metrics[f"source_fixedstep_to_refined_{label}"] = _max_gap(residual, 0)
        metrics[f"source_squared_gap_reduction_{label}"] = squared_reduction
        assert squared_reduction > (.99 if name == "float" else .9999)
    heave_limit, speed_limit = ((.008, .02) if name == "float" else (.0001, .00002))
    assert metrics["source_fixedstep_to_refined_heave"] < heave_limit
    assert metrics["source_fixedstep_to_refined_speed"] < speed_limit
    np.testing.assert_allclose(
        np.abs(fixed[:, 5] - published[:, 5]),
        np.abs(fixed[:, column] - published[:, column]),
        rtol=0, atol=1e-8,
    )

    # Python advances both bodies from the HDF5 and initial joint extension;
    # no saved source state, force, or command is passed to its solver.
    python = _physical_run(name)
    np.testing.assert_allclose(python.time[::8], published[:, 0], rtol=0, atol=1e-10)
    active = python.bodies[name]
    for label, actual, index in (
        ("heave", active.position[::8, 2], column),
        ("speed", active.velocity[::8, 2], column + 1),
    ):
        fixed_gap = _max_gap(actual, fixed[:, index])
        published_gap = _max_gap(actual, published[:, index])
        metrics[f"python_to_fixedstep_{label}"] = fixed_gap
        metrics[f"python_to_published_{label}"] = published_gap
        limit = ((.006 if label == "heave" else .02) if name == "float"
                 else (.003 if label == "heave" else .0006))
        assert fixed_gap < limit
        assert fixed_gap < .2 * published_gap
    assert python.wave_elevation is None
    print(f"RM3 PTO extension {name} fixed-step audit: {json.dumps(metrics, sort_keys=True)}")

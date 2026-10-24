"""Audit a derived, heave-only source added-mass loop break over 150 s."""

import json
import os
from pathlib import Path

import numpy as np
import pytest

from wecsim import RegularCICWave, RegularWave, WEC, WorldPoint


ROOT = Path(__file__).resolve().parents[1]
APPLICATIONS = os.environ.get("WEC_SIM_APPLICATIONS_DIR")
RECORDS = os.environ.get("WEC_SIM_ELLIPSOID_MASS_CORRECTION_DIR")
SELECTED_CASE = os.environ.get("WEC_SIM_ELLIPSOID_CASE")
ORIGINAL = (ROOT / "tests/fixtures/matlab/ellipsoid_ode45_nonlinear_clock"
            / "source_maxstep_common_grid.npz")
REFINED = (ROOT / "tests/fixtures/matlab/ellipsoid_ode45_maxstep_only"
           / "refined_005_common_grid.npz")
pytestmark = pytest.mark.skipif(
    not (APPLICATIONS and RECORDS),
    reason="fresh pinned derived MATLAB mass-correction run absent",
)

# Reuse the established refined-step physical gates, with the original
# published 0.05 s source gate still exercised separately and unchanged.
PHYSICAL_LIMITS = {
    "Regular": (.003, .004, 5000, 4000),
    "RegularCIC": (.0035, .0045, 6000, 5000),
}


def _read(folder, mode):
    shapes = {"body": (3001, 55), "pto": (3001, 25),
              "applied_added": (3001, 7),
              "wave": (3001, 2), "mass": (1, 7)}
    result = {}
    for name, shape in shapes.items():
        value = np.loadtxt(folder / f"{mode}_{name}.csv",
                           delimiter=",", ndmin=2)
        assert value.shape == shape and np.isfinite(value).all(), (mode, name)
        result[name] = value
    time = np.arange(3001) * .05
    for name in ("body", "pto", "wave", "applied_added"):
        np.testing.assert_allclose(result[name][:, 0], time,
                                   rtol=0, atol=1e-10)
    factor, trace, a33, rigid, adjusted, residual, step = result["mass"][0]
    np.testing.assert_allclose(factor, a33 / trace, rtol=0, atol=1e-13)
    np.testing.assert_allclose(adjusted, rigid + a33, rtol=0, atol=1e-7)
    assert abs(residual) < 1e-7
    assert abs(step - (.05 if mode == "published_clock" else .005)) < 1e-12
    assert np.max(np.abs(result["applied_added"][:, 3])) < 1e-5
    return result


def _error(left, right):
    assert left.shape == right.shape
    return float(np.max(np.abs(left - right)))


@pytest.mark.parametrize(
    "wave_case", [SELECTED_CASE] if SELECTED_CASE else ["Regular", "RegularCIC"],
)
def test_derived_heave_mass_shift_and_original_clock_gap(wave_case):
    assert wave_case in PHYSICAL_LIMITS
    folder = Path(RECORDS) / wave_case
    published = _read(folder, "published_clock")
    refined = _read(folder, "refined_clock")
    key = "regular" if wave_case == "Regular" else "regular_cic"
    with (np.load(ORIGINAL, allow_pickle=False) as original,
          np.load(REFINED, allow_pickle=False) as original_refined):
        unchanged_heave = original[f"{key}_published_heave"]
        unchanged_speed = original[f"{key}_published_speed"]
        unchanged_wave = original[f"{key}_published_wave"]
        fine_heave = original_refined[f"{key}_heave"]
        fine_speed = original_refined[f"{key}_speed"]
        fine_wave = original_refined[f"{key}_wave"]
    for record in (published, refined):
        assert _error(record["wave"][:, 1], unchanged_wave) < 1e-10
        assert _error(record["wave"][:, 1], fine_wave) < 1e-10

    app = Path(APPLICATIONS) / "Nonlinear_Hydro"
    wec = WEC("derived ellipsoid heave mass correction")
    body = wec.body(
        "ellipsoid", app / "hydroData/ellipsoid.h5",
        mass="equilibrium", inertia=(1.375264e6, 1.375264e6, 1.341721e6),
        geometry_file=app / "geometry/elipsoid.stl",
        nonlinear_hydro="instantaneous", drag_coefficient=1,
        drag_area=np.pi * 25,
    )
    wec.coordinate("heave", body.move("heave"))
    wec.pto("PTO1", WorldPoint(0, 0, -12.5), body.at(0, 0, 0),
            damping=1_200_000)
    wave = RegularWave(4, 6) if wave_case == "Regular" else RegularCICWave(4, 6)
    python = wec.run(
        wave, dt=.05, end_time=150, ramp_time=50, rho=1025,
        radiation_memory=60 if wave_case == "RegularCIC" else None,
    )
    position = python.bodies["ellipsoid"].position[:, 2]
    speed = python.bodies["ellipsoid"].velocity[:, 2]
    pto_force = python.ptos["PTO1"].force
    power = python.ptos["PTO1"].absorbed_power
    source = refined["body"]
    pto = refined["pto"]
    physical_gaps = (
        _error(position, source[:, 3]),
        _error(speed, source[:, 9]),
        _error(pto_force, pto[:, 15]),
        _error(power, -pto[:, 15] * pto[:, 9]),
    )
    metrics = {
        "case": wave_case,
        "corrected_mass_factor": published["mass"][0, 0],
        "corrected_published_to_refined_heave_m": _error(
            published["body"][:, 3], refined["body"][:, 3]),
        "unchanged_published_to_refined_heave_m": _error(
            unchanged_heave, fine_heave),
        "corrected_published_vs_unchanged_heave_m": _error(
            published["body"][:, 3], unchanged_heave),
        "corrected_refined_vs_unchanged_heave_m": _error(
            refined["body"][:, 3], fine_heave),
        "corrected_published_vs_unchanged_speed_mps": _error(
            published["body"][:, 9], unchanged_speed),
        "corrected_refined_vs_unchanged_speed_mps": _error(
            refined["body"][:, 9], fine_speed),
        "python_vs_corrected_refined_heave_m": physical_gaps[0],
        "python_vs_corrected_refined_speed_mps": physical_gaps[1],
        "python_vs_corrected_refined_pto_force_n": physical_gaps[2],
        "python_vs_corrected_refined_power_w": physical_gaps[3],
        "corrected_nonheave_speed_max_mps": max(
            float(np.max(np.abs(record["body"][:, [7, 8, 10, 11, 12]])))
            for record in (published, refined)
        ),
    }
    print("derived ellipsoid mass correction:", json.dumps(metrics, sort_keys=True))
    for actual, limit, name in zip(
        physical_gaps, PHYSICAL_LIMITS[wave_case],
        ("heave position", "heave speed", "PTO force", "absorbed power"),
    ):
        assert actual < limit, f"corrected refined {name}: {actual} >= {limit}"

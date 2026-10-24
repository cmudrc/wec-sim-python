"""Portable 150 s gate for a derived, heave-only source mass-loop correction."""

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
ARCHIVE = ROOT / "tests/fixtures/matlab/ellipsoid_ode45_mass_correction"
ORIGINAL = (ROOT / "tests/fixtures/matlab/ellipsoid_ode45_nonlinear_clock"
            / "source_maxstep_common_grid.npz")
REFINED = (ROOT / "tests/fixtures/matlab/ellipsoid_ode45_maxstep_only"
           / "refined_005_common_grid.npz")

# The 0.005 s limits are the established independent refined-step gates.
REFINED_PHYSICAL_LIMITS = {
    "Regular": (.003, .004, 5000, 4000),
    "RegularCIC": (.0035, .0045, 6000, 5000),
}
# The unchanged published-input test treats these as ode45 solver-envelope
# limits, not as a claim of physical parity with ode45 internals.
PUBLISHED_SOLVER_ENVELOPE_LIMITS = {
    "Regular": (.022, .032, 38_000, 30_000),
    "RegularCIC": (.023, .034, 40_000, 33_000),
}
# Established source step-stability screens from the prior refinement gate.
SOURCE_STEP_LIMITS = {
    "Regular": (.003, .00325, 4000, 4000),
    "RegularCIC": (.00375, .00425, 5000, 5000),
}


class CorrectedPublishedPhysicalGap(AssertionError):
    """The corrected 0.05 s source still misses physical parity limits."""


class CorrectedSourceStepGap(AssertionError):
    """The corrected source still changes beyond its step-stability screen."""


def _csv(folder, mode, name, shape):
    values = np.loadtxt(folder / f"{mode}_{name}.csv",
                        delimiter=",", ndmin=2)
    assert values.shape == shape and np.isfinite(values).all(), (mode, name)
    return values


def _read(wave_case, mode):
    if RECORDS:
        folder = Path(RECORDS) / wave_case
        body = _csv(folder, mode, "body", (3001, 55))
        pto = _csv(folder, mode, "pto", (3001, 25))
        wave = _csv(folder, mode, "wave", (3001, 2))
        applied = _csv(folder, mode, "applied_added", (3001, 7))
        mass = _csv(folder, mode, "mass", (1, 12))[0]
        record = {
            "time": body[:, 0], "heave": body[:, 3],
            "speed": body[:, 9], "pto_force": pto[:, 15],
            "pto_speed": pto[:, 9], "wave": wave[:, 1],
            "applied_heave": applied[:, 3],
            "inactive_speed": body[:, [7, 8, 10, 11, 12]],
            "inactive_acceleration": body[:, [49, 50, 52, 53, 54]],
            "mass": mass,
        }
        for name, values in (("pto", pto), ("wave", wave),
                             ("applied", applied)):
            np.testing.assert_allclose(values[:, 0], record["time"],
                                       rtol=0, atol=1e-10, err_msg=name)
    else:
        filename = "regular.npz" if wave_case == "Regular" else "regular_cic.npz"
        with np.load(ARCHIVE / filename, allow_pickle=False) as source:
            assert str(source["case"]) == wave_case
            record = {name: source[f"{mode}_{name}"] for name in (
                "time", "heave", "speed", "pto_force", "pto_speed",
                "wave", "applied_heave", "inactive_speed",
                "inactive_acceleration", "mass",
            )}
    time = np.arange(3001) * .05
    np.testing.assert_allclose(record["time"], time, rtol=0, atol=1e-10)
    for name in ("heave", "speed", "pto_force", "pto_speed", "wave",
                 "applied_heave"):
        assert record[name].shape == (3001,) and np.isfinite(record[name]).all()
    for name in ("inactive_speed", "inactive_acceleration"):
        assert record[name].shape == (3001, 5)
        assert np.isfinite(record[name]).all()
    factor, trace, a33, rigid, adjusted, residual, step = record["mass"][:7]
    assert record["mass"].shape == (12,)
    np.testing.assert_allclose(factor, a33 / trace, rtol=0, atol=1e-13)
    np.testing.assert_allclose(adjusted, rigid + a33, rtol=0, atol=1e-7)
    assert abs(residual) < 1e-7
    assert abs(step - (.05 if mode == "published_clock" else .005)) < 1e-12
    assert np.max(np.abs(record["inactive_speed"])) < 1e-5
    assert np.max(np.abs(record["inactive_acceleration"])) < 1e-5
    cross_force = record["inactive_acceleration"] @ record["mass"][7:]
    assert np.max(np.abs(cross_force)) < 1e-5
    assert np.max(np.abs(record["applied_heave"])) < 1e-5
    return record


def _error(left, right):
    assert left.shape == right.shape
    return float(np.max(np.abs(left - right)))


def _source_power(record):
    return -record["pto_force"] * record["pto_speed"]


def _gaps(position, speed, force, power, source):
    return (
        _error(position, source["heave"]),
        _error(speed, source["speed"]),
        _error(force, source["pto_force"]),
        _error(power, _source_power(source)),
    )


def _check_gaps(gaps, limits, error_type, label):
    names = ("heave position", "heave speed", "PTO force", "absorbed power")
    failures = [f"{name}: {actual:.8g} >= {limit:.8g}"
                for name, actual, limit in zip(names, gaps, limits)
                if not np.isfinite(actual) or actual >= limit]
    if failures:
        raise error_type(f"{label}: " + "; ".join(failures))


@pytest.fixture(scope="module", params=(
    [SELECTED_CASE] if SELECTED_CASE else ["Regular", "RegularCIC"]
))
def paired_case(request):
    wave_case = request.param
    assert wave_case in REFINED_PHYSICAL_LIMITS
    published = _read(wave_case, "published_clock")
    refined = _read(wave_case, "refined_clock")
    key = "regular" if wave_case == "Regular" else "regular_cic"
    with (np.load(ORIGINAL, allow_pickle=False) as original,
          np.load(REFINED, allow_pickle=False) as original_refined):
        unchanged_heave = original[f"{key}_published_heave"]
        unchanged_speed = original[f"{key}_published_speed"]
        unchanged_wave = original[f"{key}_published_wave"]
        fine_heave = original_refined[f"{key}_heave"]
        fine_speed = original_refined[f"{key}_speed"]
        fine_wave = original_refined[f"{key}_wave"]
    for source in (published, refined):
        assert _error(source["wave"], unchanged_wave) < 1e-10
        assert _error(source["wave"], fine_wave) < 1e-10
    if APPLICATIONS:
        app = Path(APPLICATIONS) / "Nonlinear_Hydro"
        hydro = app / "hydroData/ellipsoid.h5"
        geometry = app / "geometry/elipsoid.stl"
    else:
        hydro = ARCHIVE / "ellipsoid.h5"
        geometry = (ROOT / "tests/test_objects/test_bodyclass/testData"
                    / "body_10_test/elipsoid.stl")
    wec = WEC("derived ellipsoid heave mass correction")
    body = wec.body(
        "ellipsoid", hydro,
        mass="equilibrium", inertia=(1.375264e6, 1.375264e6, 1.341721e6),
        geometry_file=geometry, nonlinear_hydro="instantaneous",
        drag_coefficient=1, drag_area=np.pi * 25,
    )
    wec.coordinate("heave", body.move("heave"))
    wec.pto("PTO1", WorldPoint(0, 0, -12.5), body.at(0, 0, 0),
            damping=1_200_000)
    wave = RegularWave(4, 6) if wave_case == "Regular" else RegularCICWave(4, 6)
    python = wec.run(
        wave, dt=.05, end_time=150, ramp_time=50, rho=1025,
        radiation_memory=60 if wave_case == "RegularCIC" else None,
    )
    return {
        "case": wave_case, "published": published, "refined": refined,
        "unchanged_heave": unchanged_heave, "unchanged_speed": unchanged_speed,
        "fine_heave": fine_heave, "fine_speed": fine_speed,
        "position": python.bodies["ellipsoid"].position[:, 2],
        "speed": python.bodies["ellipsoid"].velocity[:, 2],
        "force": python.ptos["PTO1"].force,
        "power": python.ptos["PTO1"].absorbed_power,
    }


def test_derived_heave_mass_shift_and_refined_physical_gate(paired_case):
    d = paired_case
    case, published, refined = d["case"], d["published"], d["refined"]
    refined_gaps = _gaps(d["position"], d["speed"], d["force"],
                         d["power"], refined)
    published_gaps = _gaps(d["position"], d["speed"], d["force"],
                           d["power"], published)
    published_solver_envelope_limits = PUBLISHED_SOLVER_ENVELOPE_LIMITS[case]
    source_step_gaps = _gaps(published["heave"], published["speed"],
                             published["pto_force"], _source_power(published),
                             refined)
    metrics = {
        "case": case,
        "corrected_mass_factor": float(published["mass"][0]),
        "corrected_published_to_refined_heave_m": source_step_gaps[0],
        "corrected_published_to_refined_speed_mps": source_step_gaps[1],
        "corrected_published_to_refined_pto_force_n": source_step_gaps[2],
        "corrected_published_to_refined_power_w": source_step_gaps[3],
        "unchanged_published_to_refined_heave_m": _error(
            d["unchanged_heave"], d["fine_heave"]),
        "corrected_published_vs_unchanged_heave_m": _error(
            published["heave"], d["unchanged_heave"]),
        "corrected_refined_vs_unchanged_heave_m": _error(
            refined["heave"], d["fine_heave"]),
        "corrected_published_vs_unchanged_speed_mps": _error(
            published["speed"], d["unchanged_speed"]),
        "corrected_refined_vs_unchanged_speed_mps": _error(
            refined["speed"], d["fine_speed"]),
        "python_vs_corrected_refined_heave_m": refined_gaps[0],
        "python_vs_corrected_refined_speed_mps": refined_gaps[1],
        "python_vs_corrected_refined_pto_force_n": refined_gaps[2],
        "python_vs_corrected_refined_power_w": refined_gaps[3],
        "python_vs_corrected_published_heave_m": published_gaps[0],
        "python_vs_corrected_published_speed_mps": published_gaps[1],
        "python_vs_corrected_published_pto_force_n": published_gaps[2],
        "python_vs_corrected_published_power_w": published_gaps[3],
        "corrected_published_exceeds_original_solver_envelope": any(
            gap >= limit for gap, limit in zip(
                published_gaps, published_solver_envelope_limits)),
        "corrected_nonheave_speed_max_mps": max(
            float(np.max(np.abs(source["inactive_speed"])))
            for source in (published, refined)),
        "corrected_nonheave_acceleration_max": max(
            float(np.max(np.abs(source["inactive_acceleration"])))
            for source in (published, refined)),
    }
    print("derived ellipsoid mass correction:", json.dumps(metrics, sort_keys=True))
    _check_gaps(refined_gaps, REFINED_PHYSICAL_LIMITS[case], AssertionError,
                "derived refined source against independent Python")


@pytest.mark.xfail(
    strict=True, raises=CorrectedPublishedPhysicalGap,
    reason="derived 0.05 s source still misses refined physical limits",
)
def test_derived_corrected_published_clock_physical_gate(paired_case):
    d = paired_case
    gaps = _gaps(d["position"], d["speed"], d["force"], d["power"],
                  d["published"])
    _check_gaps(gaps, REFINED_PHYSICAL_LIMITS[d["case"]],
                CorrectedPublishedPhysicalGap, "derived published clock")


@pytest.mark.xfail(
    strict=True, raises=CorrectedSourceStepGap,
    reason="derived 0.05 s source has not passed its step-stability screen",
)
def test_derived_corrected_source_step_stability(paired_case):
    d = paired_case
    published, refined = d["published"], d["refined"]
    gaps = _gaps(published["heave"], published["speed"],
                  published["pto_force"], _source_power(published), refined)
    _check_gaps(gaps, SOURCE_STEP_LIMITS[d["case"]],
                CorrectedSourceStepGap, "derived source 0.05 to 0.005 s")

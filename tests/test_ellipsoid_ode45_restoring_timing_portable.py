"""Portable 150 s screen for the copied ellipsoid source force-timing arms.

The unchanged published 0.05 s source remains the separate original-input
gate. These edited-source arms are diagnostics; their inherited physical
motion, PTO, and work limits still fail.
"""

from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import trapezoid

from test_ellipsoid_ode45_convergence import PYTHON_LIMITS, _python
from wecsim.hydrodynamics.nonlinearHydro import HeaveMeshHydro


ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests/fixtures/matlab/ellipsoid_ode45_restoring_timing"
SOURCE_INPUT = ROOT / "tests/fixtures/matlab/ellipsoid_ode45_convergence"
ORIGINAL = (ROOT / "tests/fixtures/matlab/ellipsoid_ode45_nonlinear_clock"
            / "source_maxstep_common_grid.npz")
MODES = ("baseline", "continuous_restoring", "continuous_nonlinear")
CHANNELS = ("heave", "speed", "pto_stroke", "pto_speed", "pto_force",
            "pto_power")
FIELDS = ("time", "heave", "speed", "restoring", "excitation",
          "acceleration", "pto_stroke", "pto_speed", "pto_force",
          "pto_power", "wave", "settings")
SAMPLE_PERIODS = {
    "baseline": (.05, .05, .05),
    "continuous_restoring": (0., .05, .05),
    "continuous_nonlinear": (0., 0., 0.),
}


def _peak_gap(actual, expected):
    assert actual.shape == expected.shape
    assert np.isfinite(actual).all() and np.isfinite(expected).all()
    return float(np.max(np.abs(actual - expected)))


@pytest.mark.parametrize("case", ("Regular", "RegularCIC"))
def test_full_duration_copied_source_force_timing_diagnostic(case):
    name = "regular_cic" if case == "RegularCIC" else "regular"
    with np.load(FIXTURE / f"{name}.npz", allow_pickle=False) as archive:
        expected_keys = {"case"} | {
            f"{mode}_{field}" for mode in MODES for field in FIELDS
        }
        assert set(archive.files) == expected_keys
        assert str(archive["case"]) == case
        records = {
            mode: {field: archive[f"{mode}_{field}"] for field in FIELDS}
            for mode in MODES
        }

    times = np.arange(3001) * .05
    for mode, record in records.items():
        np.testing.assert_allclose(record["settings"],
                                   [.05, .05, .05, .05, .05, 150,
                                    *SAMPLE_PERIODS[mode]],
                                   rtol=0, atol=1e-12)
        assert record["settings"].shape == (9,)
        for field in FIELDS[:-1]:
            assert record[field].shape == (3001,), (case, mode, field)
            assert np.isfinite(record[field]).all(), (case, mode, field)
        np.testing.assert_allclose(record["time"], times, rtol=0, atol=1e-9)
        assert _peak_gap(record["pto_speed"], record["speed"]) < 1e-9
        assert _peak_gap(record["pto_stroke"], record["heave"] + 2) < 1e-10
        assert _peak_gap(record["pto_force"],
                         -1_200_000 * record["pto_speed"]) < 1e-6
        assert _peak_gap(record["pto_power"],
                         -record["pto_force"] * record["pto_speed"]) < 1e-6
        assert _peak_gap(record["wave"], records["baseline"]["wave"]) < 1e-10

    # The fresh copied-model baseline reproduces the unchanged published run,
    # independently archived before this timing intervention.
    with np.load(ORIGINAL, allow_pickle=False) as original:
        baseline = records["baseline"]
        for field, limit in (("time", 1e-10), ("heave", 1e-7),
                             ("speed", 1e-7), ("restoring", .1),
                             ("excitation", .1), ("pto_force", .1),
                             ("wave", 1e-10)):
            assert _peak_gap(
                baseline[field], original[f"{name}_published_{field}"
                                          ]) < limit, (case, field)
        np.testing.assert_allclose(baseline["settings"][:5],
                                   original[f"{name}_published_settings"],
                                   rtol=0, atol=1e-12)

    mesh = HeaveMeshHydro.from_stl(
        SOURCE_INPUT / "geometry/elipsoid.stl", center_z=-2, rho=1025,
        gravity=9.81, depth=70, period=6, height=4, ramp_time=50,
        mass=None, drag_coefficient=1, drag_area=np.pi * 25,
    )
    for mode, record in records.items():
        mesh_force = np.array([
            mesh.forces(t, z + 2, v)[0]
            for t, z, v in zip(record["time"], record["heave"],
                               record["speed"], strict=True)
        ])
        applied = -record["restoring"]
        previous_error = _peak_gap(applied[1:], mesh_force[:-1])
        current_error = _peak_gap(applied[1:], mesh_force[1:])
        if mode == "baseline":
            assert previous_error < 1e-6 and current_error > 500
        else:
            assert current_error < 1e-6 and previous_error > 500

    # Python advances its own dynamics using the archived source-generated
    # HDF5/STL and published wave/PTO input; no source motion or force is fed in.
    python = _python(case, .05, SOURCE_INPUT)
    python_body = python.bodies["ellipsoid"]
    python_pto = python.ptos["PTO1"]
    physical = {
        "heave": python_body.position[:, 2],
        "speed": python_body.velocity[:, 2],
        "pto_stroke": python_pto.stroke,
        "pto_speed": python_pto.velocity,
        "pto_force": python_pto.force,
        "pto_power": python_pto.absorbed_power,
    }
    assert _peak_gap(python.time, times) < 1e-9
    assert _peak_gap(python.wave_elevation, baseline["wave"]) < 1e-10
    python_work = float(trapezoid(physical["pto_power"], python.time))
    assert python_work > 1e3

    gaps = {
        mode: {field: _peak_gap(record[field], physical[field])
               for field in CHANNELS}
        for mode, record in records.items()
    }
    work_gaps = {
        mode: abs(float(trapezoid(record["pto_power"], times))
                  / python_work - 1)
        for mode, record in records.items()
    }
    pointwise_limits = PYTHON_LIMITS[case]
    for mode in MODES[1:]:
        # These are regression screens for the observed source improvement,
        # not prospective acceptance limits. Both edited source routes still
        # fail every inherited full-duration physical gate.
        for field in CHANNELS:
            assert gaps[mode][field] < .5 * gaps["baseline"][field], (
                case, mode, field, gaps[mode][field], gaps["baseline"][field]
            )
        assert work_gaps[mode] < .5 * work_gaps["baseline"]
        for field, limit in zip(
            ("heave", "speed", "pto_force", "pto_power"),
            pointwise_limits, strict=True,
        ):
            assert gaps[mode][field] > limit, (case, mode, field)
        assert gaps[mode]["pto_stroke"] > pointwise_limits[0]
        assert gaps[mode]["pto_speed"] > pointwise_limits[1]
        assert work_gaps[mode] > .001
    print(f"{case} edited-source peak gaps: {gaps}; work fractions: {work_gaps}")

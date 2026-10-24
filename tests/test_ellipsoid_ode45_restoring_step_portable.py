"""Portable copied-model MaxStep evidence after the failed 0.0125 s screen.

The 0.00625 and 0.003125 s arms are later diagnostics. Their comparison with
Python does not qualify the unchanged published 0.05 s source input.
"""

from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import trapezoid

from test_ellipsoid_ode45_convergence import _python
from tools.check_ellipsoid_ode45_restoring_step import (
    PHYSICAL_LIMITS, SOURCE_LIMITS,
)
from wecsim.hydrodynamics.nonlinearHydro import HeaveMeshHydro


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "tests/fixtures/matlab/ellipsoid_ode45_restoring_step"
BASE = ROOT / "tests/fixtures/matlab/ellipsoid_ode45_restoring_timing"
SOURCE_INPUT = ROOT / "tests/fixtures/matlab/ellipsoid_ode45_convergence"
STEPS = (("0p025", .025), ("0p0125", .0125),
         ("0p00625", .00625), ("0p003125", .003125))
FIELDS = ("heave", "speed", "restoring", "pto_stroke", "pto_speed",
          "pto_force", "pto_power", "wave", "settings")
CHANNELS = ("heave", "speed", "pto_force", "pto_power")
PREDECLARED_FAILURES = {
    "Regular": {"speed", "pto_force", "pto_power", "signed_work"},
    "RegularCIC": {"heave", "speed", "pto_force", "pto_power",
                   "signed_work"},
}


def _peak(left, right):
    assert left.shape == right.shape
    assert np.isfinite(left).all() and np.isfinite(right).all()
    return float(np.max(np.abs(left - right)))


@pytest.mark.parametrize("case", ("Regular", "RegularCIC"))
def test_copied_model_maxstep_continuation(case):
    name = "regular_cic" if case == "RegularCIC" else "regular"
    with np.load(BASE / f"{name}.npz", allow_pickle=False) as archive:
        assert str(archive["case"]) == case
        base = {field: archive[f"continuous_nonlinear_{field}"]
                for field in FIELDS}
        base["time"] = archive["continuous_nonlinear_time"]
    with np.load(ARCHIVE / f"{name}.npz", allow_pickle=False) as archive:
        expected = {"case", "time"} | {
            f"maxstep_{tag}_{field}" for tag, _ in STEPS for field in FIELDS
        }
        assert set(archive.files) == expected
        assert str(archive["case"]) == case
        time = archive["time"]
        records = {
            tag: {field: archive[f"maxstep_{tag}_{field}"]
                  for field in FIELDS}
            for tag, _ in STEPS
        }

    assert time.shape == (3001,)
    np.testing.assert_allclose(time, np.arange(3001) * .05,
                               rtol=0, atol=1e-9)
    np.testing.assert_allclose(base["time"], time, rtol=0, atol=1e-9)
    np.testing.assert_allclose(base["settings"],
                               [.05, .05, .05, .05, .05, 150, 0, 0, 0],
                               rtol=0, atol=1e-12)

    mesh = HeaveMeshHydro.from_stl(
        SOURCE_INPUT / "geometry/elipsoid.stl", center_z=-2, rho=1025,
        gravity=9.81, depth=70, period=6, height=4, ramp_time=50,
        mass=None, drag_coefficient=1, drag_area=np.pi * 25,
    )
    for tag, step in STEPS:
        record = records[tag]
        np.testing.assert_allclose(
            record["settings"], [.05, .05, .05, .05, step, 150, 0, 0, 0],
            rtol=0, atol=1e-12,
        )
        assert record["settings"].shape == (9,)
        for field in FIELDS[:-1]:
            assert record[field].shape == (3001,), (case, tag, field)
            assert np.isfinite(record[field]).all(), (case, tag, field)
        assert _peak(record["wave"], base["wave"]) < 1e-10
        assert _peak(record["pto_stroke"], record["heave"] + 2) < 1e-10
        assert _peak(record["pto_speed"], record["speed"]) < 1e-9
        assert _peak(record["pto_force"],
                     -1_200_000 * record["pto_speed"]) < 1e-6
        assert _peak(record["pto_power"],
                     -record["pto_force"] * record["pto_speed"]) < 1e-6
        expected_restoring = np.array([
            mesh.forces(t, z + 2, v)[0]
            for t, z, v in zip(time, record["heave"], record["speed"],
                               strict=True)
        ])
        assert _peak(-record["restoring"], expected_restoring) < 1e-6

    histories = [base] + [records[tag] for tag, _ in STEPS]
    source_changes = {
        field: [_peak(a[field], b[field])
                for a, b in zip(histories, histories[1:])]
        for field in CHANNELS
    }
    work = [float(trapezoid(record["pto_power"], time))
            for record in histories]
    work_changes = [abs(later - earlier) / later
                    for earlier, later in zip(work, work[1:])]

    # The .025 -> .0125 result was the prospective screen. Preserve its
    # observed failure, including the signed-work limit, as historical fact.
    failures = set()
    for field, limit in zip(CHANNELS, SOURCE_LIMITS[case], strict=True):
        coarse, screened = source_changes[field][:2]
        assert screened < .7 * coarse
        if screened >= limit:
            failures.add(field)
    assert work_changes[1] < .7 * work_changes[0]
    if work_changes[1] >= .001:
        failures.add("signed_work")
    assert failures == PREDECLARED_FAILURES[case]
    print(f"{case} predeclared 0.0125 s screen failures: {sorted(failures)}")
    print(f"{case} later source changes: {source_changes}; "
          f"signed-work fractions: {work_changes}")

    # Independent Python dynamics use only the archived source HDF5/STL and
    # configured wave/PTO. The last arm is a post-screen physical diagnostic.
    python = _python(case, .05, SOURCE_INPUT)
    body = python.bodies["ellipsoid"]
    pto = python.ptos["PTO1"]
    physical = {
        "heave": body.position[:, 2],
        "speed": body.velocity[:, 2],
        "pto_stroke": pto.stroke,
        "pto_speed": pto.velocity,
        "pto_force": pto.force,
        "pto_power": pto.absorbed_power,
    }
    assert _peak(python.time, time) < 1e-9
    assert _peak(python.wave_elevation, base["wave"]) < 1e-10
    for tag, _ in STEPS:
        assert _peak(records[tag]["wave"], python.wave_elevation) < 1e-10
    final = records["0p003125"]
    gaps = {field: _peak(final[field], values)
            for field, values in physical.items()}
    limits = dict(zip(CHANNELS, PHYSICAL_LIMITS[case], strict=True))
    limits["pto_stroke"] = limits["heave"]
    limits["pto_speed"] = limits["speed"]
    for field, limit in limits.items():
        assert gaps[field] < limit, (case, field, gaps[field], limit)
    python_work = float(trapezoid(physical["pto_power"], time))
    final_work_error = abs(work[-1] / python_work - 1)
    assert final_work_error < .001
    print(f"{case} post-screen 0.003125 s Python gaps: {gaps}; "
          f"signed-work fraction: {final_work_error}")

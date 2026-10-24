"""Portable post-screen MaxStep diagnosis without active heave delay.

These copied-model controls follow a failed prospective source-step screen.
Their refined behavior cannot qualify either unchanged published input.
"""

from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import trapezoid

from test_ellipsoid_ode45_convergence import PYTHON_LIMITS, _python
from wecsim.hydrodynamics.nonlinearHydro import HeaveMeshHydro


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = (ROOT / "tests/fixtures/matlab"
           / "ellipsoid_ode45_continuous_nonlinear_mass_control_step")
BASE = (ROOT / "tests/fixtures/matlab"
        / "ellipsoid_ode45_continuous_nonlinear_mass_control")
DELAYED = ROOT / "tests/fixtures/matlab/ellipsoid_ode45_restoring_step"
DELAYED_BASE = ROOT / "tests/fixtures/matlab/ellipsoid_ode45_restoring_timing"
SOURCE_INPUT = ROOT / "tests/fixtures/matlab/ellipsoid_ode45_convergence"
STEPS = (("0p0125", .0125), ("0p00625", .00625))
TIME = np.arange(3001) * .05
CHANNELS = ("heave", "speed", "pto_stroke", "pto_speed", "pto_force",
            "pto_power")
FIELDS = ("heave", "speed", "total", "excitation", "radiation",
          "reported_added_mass", "restoring", "viscous", "acceleration",
          "applied_heave_added", "inactive_speed", "inactive_acceleration",
          "pto_stroke", "pto_speed", "pto_force", "pto_power", "wave",
          "settings", "mass")


def _peak(left, right):
    assert left.shape == right.shape
    assert np.isfinite(left).all() and np.isfinite(right).all()
    return float(np.max(np.abs(left - right)))


@pytest.mark.parametrize("case", ("Regular", "RegularCIC"))
def test_post_screen_no_active_delay_gaps_contract_but_work_remains_outside_gate(case):
    filename = "regular_cic.npz" if case == "RegularCIC" else "regular.npz"
    with (np.load(BASE / filename, allow_pickle=False) as a,
          np.load(ARCHIVE / filename, allow_pickle=False) as fine,
          np.load(DELAYED_BASE / filename, allow_pickle=False) as b,
          np.load(DELAYED / filename, allow_pickle=False) as delayed):
        assert set(fine.files) == ({"case", "time"} | {
            f"maxstep_{tag}_{field}" for tag, _ in STEPS for field in FIELDS
        })
        assert str(fine["case"]) == str(a["case"]) == case
        np.testing.assert_allclose(fine["time"], TIME, rtol=0, atol=1e-9)
        base = {field: a[field] for field in FIELDS}
        base["time"] = a["time"]
        delayed_base = {
            field: b[f"continuous_nonlinear_{field}"]
            for field in CHANNELS + ("wave",)
        }
        controls = {
            tag: {field: fine[f"maxstep_{tag}_{field}"] for field in FIELDS}
            for tag, _ in STEPS
        }
        paired = {
            tag: {
                field: delayed[f"maxstep_{tag}_{field}"]
                for field in CHANNELS + ("wave",)
            }
            for tag, _ in STEPS
        }

    np.testing.assert_allclose(base["time"], TIME, rtol=0, atol=1e-9)
    np.testing.assert_allclose(base["settings"],
                               [.05, .05, .05, .05, .05, 150, 0, 0, 0],
                               rtol=0, atol=1e-12)
    mesh = HeaveMeshHydro.from_stl(
        SOURCE_INPUT / "geometry/elipsoid.stl", center_z=-2, rho=1025,
        gravity=9.81, depth=70, period=6, height=4, ramp_time=50,
        mass=None, drag_coefficient=1, drag_area=np.pi * 25,
    )
    assert abs(mesh.mass - base["mass"][3]) < 1e-7
    records = {"0p05": base, **controls}
    matches = {"0p05": delayed_base, **paired}
    for tag, step in STEPS:
        record = records[tag]
        np.testing.assert_allclose(record["settings"],
                                   [.05, .05, .05, .05, step, 150, 0, 0, 0],
                                   rtol=0, atol=1e-12)
        for field in FIELDS:
            shape = (3001, 5) if field.startswith("inactive_") else (
                (9,) if field == "settings" else (12,) if field == "mass"
                else (3001,)
            )
            assert record[field].shape == shape, (case, tag, field)
            assert np.isfinite(record[field]).all(), (case, tag, field)
        mass = record["mass"]
        np.testing.assert_allclose(np.delete(mass, 6),
                                   np.delete(base["mass"], 6),
                                   rtol=0, atol=1e-7)
        assert abs(mass[5]) < 1e-7 and abs(mass[6] - step) < 1e-12
        assert np.max(np.abs(record["inactive_speed"])) < 1e-9
        assert np.max(np.abs(record["inactive_acceleration"])) < 1e-9
        assert np.max(np.abs(record["inactive_acceleration"] @ mass[7:])) < 1e-6
        assert np.max(np.abs(record["applied_heave_added"])) < 1e-6
        assert np.max(np.abs(record["reported_added_mass"])) < 1e-6
        assert _peak(record["wave"], base["wave"]) < 1e-10
        assert _peak(record["wave"], matches[tag]["wave"]) < 1e-10
        assert _peak(record["pto_stroke"], record["heave"] + 2) < 1e-9
        assert _peak(record["pto_speed"], record["speed"]) < 1e-9
        assert _peak(record["pto_force"],
                     -1_200_000 * record["pto_speed"]) < 1e-6
        assert _peak(record["pto_power"],
                     -record["pto_force"] * record["pto_speed"]) < 1e-6
        expected = np.array([
            mesh.forces(t, z + 2, v)
            for t, z, v in zip(TIME, record["heave"], record["speed"],
                               strict=True)
        ])
        assert _peak(-record["restoring"], expected[:, 0]) < 1e-6
        assert _peak(-record["viscous"], expected[:, 2]) < 1e-6

    # Python's dynamics use only the pinned HDF5/STL and wave/PTO inputs.
    python = _python(case, .05, SOURCE_INPUT)
    body, pto = python.bodies["ellipsoid"], python.ptos["PTO1"]
    physical = {
        "heave": body.position[:, 2], "speed": body.velocity[:, 2],
        "pto_stroke": pto.stroke, "pto_speed": pto.velocity,
        "pto_force": pto.force, "pto_power": pto.absorbed_power,
    }
    assert _peak(python.time, TIME) < 1e-9
    assert _peak(python.wave_elevation, base["wave"]) < 1e-10
    control_gaps = {
        tag: {field: _peak(record[field], physical[field])
              for field in CHANNELS}
        for tag, record in records.items()
    }
    paired_gaps = {
        tag: {field: _peak(records[tag][field], matches[tag][field])
              for field in CHANNELS}
        for tag in records
    }
    control_work = {
        tag: float(trapezoid(record["pto_power"], TIME))
        for tag, record in records.items()
    }
    paired_work = {
        tag: float(trapezoid(matches[tag]["pto_power"], TIME))
        for tag in records
    }
    python_work = float(trapezoid(physical["pto_power"], TIME))
    assert min(*control_work.values(), *paired_work.values(), python_work) > 1e3
    work_gaps = {tag: abs(work/python_work - 1)
                 for tag, work in control_work.items()}
    paired_work_gaps = {tag: abs(control_work[tag]/paired_work[tag] - 1)
                        for tag in records}

    order = ("0p05", "0p0125", "0p00625")
    for field in CHANNELS:
        assert all(control_gaps[later][field] < control_gaps[earlier][field]
                   for earlier, later in zip(order, order[1:])), (case, field)
        assert all(paired_gaps[later][field] < paired_gaps[earlier][field]
                   for earlier, later in zip(order, order[1:])), (case, field)
    assert all(work_gaps[later] < work_gaps[earlier]
               for earlier, later in zip(order, order[1:]))
    assert all(paired_work_gaps[later] < paired_work_gaps[earlier]
               for earlier, later in zip(order, order[1:]))

    # The edited post-screen arms satisfy pointwise inherited limits, while
    # their signed work still exceeds 0.1%. The original gate is unchanged.
    limits = dict(zip(("heave", "speed", "pto_force", "pto_power"),
                      PYTHON_LIMITS[case], strict=True))
    limits["pto_stroke"] = limits["heave"]
    limits["pto_speed"] = limits["speed"]
    for tag, _ in STEPS:
        for field, limit in limits.items():
            assert control_gaps[tag][field] < limit, (
                case, tag, field, control_gaps[tag][field], limit
            )
        assert work_gaps[tag] > .001, (case, tag, work_gaps[tag])
    print(f"{case} no-delay/Python peak gaps: {control_gaps}")
    print(f"{case} no-delay/delayed peak gaps: {paired_gaps}")
    print(f"{case} no-delay/Python signed-work fractions: {work_gaps}; "
          f"no-delay/delayed fractions: {paired_work_gaps}")

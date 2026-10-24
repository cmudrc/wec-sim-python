"""Portable 150 s diagnostic for continuous nonlinear force plus mass shift.

The copied MATLAB model keeps the published 0.05 s clocks and evaluates all
three nonlinear blocks continuously. Moving exactly A_zz into the rigid body
removes active heave acceleration from the 100 ns delay. These source edits
remain separate from the unchanged published-input physical gate.
"""

from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import trapezoid

from test_ellipsoid_ode45_convergence import PYTHON_LIMITS, _python
from wecsim.hydrodynamics.nonlinearHydro import HeaveMeshHydro


ROOT = Path(__file__).resolve().parents[1]
CONTROL = (ROOT / "tests/fixtures/matlab"
           / "ellipsoid_ode45_continuous_nonlinear_mass_control")
TIMING = ROOT / "tests/fixtures/matlab/ellipsoid_ode45_restoring_timing"
SOURCE_INPUT = ROOT / "tests/fixtures/matlab/ellipsoid_ode45_convergence"
SOURCE_FIELDS = ("time", "heave", "speed", "restoring", "excitation",
                 "acceleration", "pto_stroke", "pto_speed", "pto_force",
                 "pto_power", "wave")
PHYSICAL_FIELDS = ("heave", "speed", "pto_stroke", "pto_speed", "pto_force",
                   "pto_power")
CONTROL_FIELDS = ("time", "heave", "speed", "total", "excitation",
                  "radiation", "reported_added_mass", "restoring", "viscous",
                  "acceleration", "applied_heave_added", "pto_stroke",
                  "pto_speed", "pto_force", "pto_power", "wave")


def _peak(left, right):
    assert left.shape == right.shape
    assert np.isfinite(left).all() and np.isfinite(right).all()
    return float(np.max(np.abs(left - right)))


@pytest.mark.parametrize("case", ("Regular", "RegularCIC"))
def test_no_active_delay_control_does_not_close_physical_gap(case):
    name = "regular_cic" if case == "RegularCIC" else "regular"
    with (np.load(CONTROL / f"{name}.npz", allow_pickle=False) as packed,
          np.load(TIMING / f"{name}.npz", allow_pickle=False) as archived):
        assert set(packed.files) == ({"case", "settings", "mass",
                                     "inactive_speed", "inactive_acceleration"}
                                    | set(CONTROL_FIELDS))
        assert str(packed["case"]) == case
        control = {key: packed[key] for key in packed.files if key != "case"}
        old = {field: archived[f"continuous_nonlinear_{field}"]
               for field in SOURCE_FIELDS}

    time = np.arange(3001) * .05
    for field in CONTROL_FIELDS:
        assert control[field].shape == (3001,), field
        assert np.isfinite(control[field]).all(), field
    np.testing.assert_allclose(control["time"], time, rtol=0, atol=1e-9)
    np.testing.assert_allclose(control["settings"],
                               [.05, .05, .05, .05, .05, 150, 0, 0, 0],
                               rtol=0, atol=1e-12)
    np.testing.assert_allclose(control["time"], old["time"],
                               rtol=0, atol=1e-9)
    assert _peak(control["wave"], old["wave"]) < 1e-10

    mass = control["mass"]
    assert mass.shape == (12,)
    factor, trace, added_zz, original, adjusted, active_residual, step = mass[:7]
    np.testing.assert_allclose(factor, added_zz / trace, rtol=0, atol=1e-13)
    np.testing.assert_allclose(adjusted, original + added_zz,
                               rtol=0, atol=1e-7)
    assert abs(active_residual) < 1e-7 and abs(step - .05) < 1e-12
    for field in ("inactive_speed", "inactive_acceleration"):
        assert control[field].shape == (3001, 5)
        assert np.isfinite(control[field]).all()
        assert np.max(np.abs(control[field])) < 1e-9
    assert np.max(np.abs(control["inactive_acceleration"] @ mass[7:])) < 1e-6
    assert np.max(np.abs(control["applied_heave_added"])) < 1e-6
    assert np.max(np.abs(control["reported_added_mass"])) < 1e-6

    assert _peak(control["pto_speed"], control["speed"]) < 1e-9
    assert _peak(control["pto_stroke"], control["heave"] + 2) < 1e-9
    assert _peak(control["pto_force"],
                 -1_200_000 * control["pto_speed"]) < 1e-6
    assert _peak(control["pto_power"],
                 -control["pto_force"] * control["pto_speed"]) < 1e-6

    mesh = HeaveMeshHydro.from_stl(
        SOURCE_INPUT / "geometry/elipsoid.stl", center_z=-2, rho=1025,
        gravity=9.81, depth=70, period=6, height=4, ramp_time=50,
        mass=None, drag_coefficient=1, drag_area=np.pi * 25,
    )
    assert abs(mesh.mass - original) < 1e-7
    forces = np.array([
        mesh.forces(t, z + 2, v)
        for t, z, v in zip(control["time"], control["heave"],
                           control["speed"], strict=True)
    ])
    # Source logs these forces with the opposite resisting-force sign.
    assert _peak(-control["restoring"], forces[:, 0]) < 1e-6
    assert _peak(-control["viscous"], forces[:, 2]) < 1e-6

    # These full-grid comparisons quantify the effect of the mass shift on
    # the previously archived continuous-nonlinear source arm.
    source_change = {
        field: _peak(control[field], old[field])
        for field in SOURCE_FIELDS if field not in ("time", "wave")
    }
    assert source_change["heave"] > .003
    assert source_change["speed"] > .007
    assert source_change["restoring"] > 9000
    assert source_change["excitation"] > 8000
    assert source_change["pto_force"] > 8000

    # Python independently advances its own 150 s dynamics from the pinned
    # HDF5/STL and wave/PTO inputs; no source state or force is fed into it.
    python = _python(case, .05, SOURCE_INPUT)
    python_body, python_pto = (python.bodies["ellipsoid"],
                               python.ptos["PTO1"])
    physical = {
        "heave": python_body.position[:, 2],
        "speed": python_body.velocity[:, 2],
        "pto_stroke": python_pto.stroke,
        "pto_speed": python_pto.velocity,
        "pto_force": python_pto.force,
        "pto_power": python_pto.absorbed_power,
    }
    assert _peak(control["time"], python.time) < 1e-9
    assert _peak(control["wave"], python.wave_elevation) < 1e-10
    old_gaps = {field: _peak(old[field], physical[field])
                for field in PHYSICAL_FIELDS}
    control_gaps = {field: _peak(control[field], physical[field])
                    for field in PHYSICAL_FIELDS}
    limits = PYTHON_LIMITS[case]
    for field, limit in zip(("heave", "speed", "pto_force", "pto_power"),
                            limits, strict=True):
        assert control_gaps[field] > limit, (case, field, control_gaps[field])
    assert control_gaps["pto_stroke"] > limits[0]
    assert control_gaps["pto_speed"] > limits[1]
    for field in PHYSICAL_FIELDS:
        assert control_gaps[field] > old_gaps[field], (
            case, field, control_gaps[field], old_gaps[field]
        )

    old_work = float(trapezoid(old["pto_power"], time))
    control_work = float(trapezoid(control["pto_power"], time))
    python_work = float(trapezoid(physical["pto_power"], time))
    assert min(old_work, control_work, python_work) > 1e3
    old_work_gap = abs(old_work / python_work - 1)
    control_work_gap = abs(control_work / python_work - 1)
    assert control_work_gap > old_work_gap > .001
    print(f"{case} mass-control source change: {source_change}")
    print(f"{case} archived/Python peak gaps: {old_gaps}")
    print(f"{case} control/Python peak gaps: {control_gaps}")
    print(f"{case} source/Python signed-work fractions: "
          f"{old_work_gap}, {control_work_gap}")

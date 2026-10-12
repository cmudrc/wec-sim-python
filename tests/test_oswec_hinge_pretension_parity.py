"""Pair a configured OSWEC hinge spring and pretension with pinned WEC-Sim."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim import RegularWave, WEC, WorldPoint


APPLICATIONS = os.environ.get("WEC_SIM_APPLICATIONS_DIR")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
pytestmark = pytest.mark.skipif(
    not (APPLICATIONS and REFERENCE),
    reason="configured OSWEC MATLAB output and HDF5 not provided",
)


def _max_error(actual, expected, limit, label):
    assert actual.shape == expected.shape, label
    assert np.isfinite(actual).all() and np.isfinite(expected).all(), label
    error = float(np.max(np.abs(actual - expected)))
    print(f"{label}: {error:.6g} (gate {limit:.6g})")
    assert error < limit, label


def test_configured_hinge_spring_and_pretension_against_matlab():
    source = Path(REFERENCE)
    parameters = np.loadtxt(
        source / "OSWEC_Nonhydro_Pretension_parameters.csv",
        delimiter=",", ndmin=1,
    )
    np.testing.assert_allclose(
        parameters, [50_000, 12_000, -5_000, 0.1], rtol=0, atol=1e-12,
    )
    bodies = [np.loadtxt(
        source / f"OSWEC_Nonhydro_Pretension_body{number}.csv",
        delimiter=",",
    ) for number in (1, 2)]
    pto = np.loadtxt(
        source / "OSWEC_Nonhydro_Pretension_pto1.csv", delimiter=",",
    )
    assert all(values.shape == (4001, 25) for values in (*bodies, pto))
    assert abs(pto[0, 17]) > 1_000  # source spring is loaded at startup

    hydro = (Path(APPLICATIONS)
             / "_Common_Input_Files/OSWEC/hydroData/oswec.h5").resolve()
    wec = WEC("OSWEC configured hinge PTO")
    flap = wec.body("flap", hydro, mass=127_000,
                    inertia=(0, 1.85e6, 0))
    base = wec.fixed_body(
        "base", center_gravity=(0, 0, -10.9), mass=999,
        inertia=(1, 1, 1),
    )
    wec.fixed_hinge(
        flap, base, location=WorldPoint(0, 0, -10),
        pto_location=WorldPoint(0, 0, -8.9), damping=12_000,
        stiffness=50_000, pretension=-5_000,
    )
    result = wec.run(
        RegularWave(2.5, 8), dt=0.1, end_time=400, ramp_time=100,
    )
    assert result.case["pto"]["equilibrium_position"] == 0.1
    _max_error(result.time, pto[:, 0], 1e-10, "time")
    for name, expected in zip(("flap", "base"), bodies):
        _max_error(result.bodies[name].position,
                   expected[:, 1:7],
                   0.04 if name == "flap" else 1e-10,
                   f"{name} position")
        _max_error(result.bodies[name].velocity,
                   expected[:, 7:13],
                   0.03 if name == "flap" else 1e-10,
                   f"{name} velocity")

    history = result.ptos["hinge"]
    _max_error(history.stroke, pto[:, 5], 0.008, "PTO angle")
    _max_error(history.velocity, pto[:, 11], 0.008, "PTO angular speed")
    _max_error(history.force, pto[:, 17], 1_000, "PTO torque")
    _max_error(history.force,
               -12_000 * history.velocity
               - 50_000 * (history.stroke - 0.1),
               1e-8, "Python PTO torque law")
    _max_error(pto[:, 17],
               -12_000 * pto[:, 11] - 50_000 * (pto[:, 5] - 0.1),
               1e-5, "MATLAB PTO torque law")
    _max_error(pto[:, 23], pto[:, 17] * pto[:, 11],
               1e-5, "MATLAB PTO mechanical power")
    _max_error(-history.force * history.velocity, -pto[:, 23],
               1_000, "PTO mechanical power")

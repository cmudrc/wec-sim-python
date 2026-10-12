"""Pair a configured RM3 PTO spring and pretension with pinned WEC-Sim."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim import RegularWave, WEC


HYDRO = os.environ.get("WEC_SIM_RM3_H5")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
pytestmark = pytest.mark.skipif(
    not (HYDRO and REFERENCE),
    reason="configured RM3 MATLAB output and HDF5 not provided",
)


def _max_error(actual, expected, limit, label):
    actual = np.asarray(actual)
    expected = np.asarray(expected)
    assert actual.shape == expected.shape, label
    assert np.isfinite(actual).all() and np.isfinite(expected).all(), label
    error = float(np.max(np.abs(actual - expected)))
    assert error < limit, f"{label}: {error:.6g} exceeds {limit}"


def test_rm3_configured_spring_and_pretension_against_matlab():
    source = Path(REFERENCE)
    parameters = np.loadtxt(
        source / "RM3_PTO_PRETENSION_parameters.csv", delimiter=",", ndmin=1,
    )
    np.testing.assert_allclose(parameters, [50_000, 5_000, -0.1], rtol=0, atol=1e-12)
    bodies = [
        np.loadtxt(source / f"RM3_PTO_PRETENSION_Pretension_body{number}.csv",
                   delimiter=",")
        for number in (1, 2)
    ]
    pto = np.loadtxt(
        source / "RM3_PTO_PRETENSION_Pretension_pto1.csv", delimiter=",",
    )
    assert all(values.shape == (4001, 25) for values in (*bodies, pto))
    assert abs(pto[0, 15]) > 1_000  # the source spring is loaded at startup

    wec = WEC("RM3 configured PTO")
    float_body = wec.body("float", Path(HYDRO).resolve(),
                          inertia=(0, 21_306_090.66, 0))
    spar = wec.body("spar", Path(HYDRO).resolve(),
                    inertia=(0, 94_407_091.24, 0))
    wec.floating_joint(
        float_body, spar, damping=1_200_000, stiffness=50_000,
        pretension=5_000,
    )
    result = wec.run(
        RegularWave(2.5, 8), dt=0.1, end_time=400, ramp_time=100,
    )
    assert result.case["pto"]["equilibrium_position"] == -0.1
    _max_error(result.time, pto[:, 0], 1e-10, "time")
    for name, expected, position_limits, speed_limits in (
        ("float", bodies[0], (0.035, 0.006, 0.001), (0.025, 0.005, 0.001)),
        ("spar", bodies[1], (0.045, 0.0025, 0.001), (0.035, 0.002, 0.001)),
    ):
        for dof, position_limit, speed_limit in zip(
            (0, 2, 4), position_limits, speed_limits,
        ):
            _max_error(result.bodies[name].position[:, dof],
                       expected[:, 1 + dof], position_limit,
                       f"{name} position DOF {dof}")
            _max_error(result.bodies[name].velocity[:, dof],
                       expected[:, 7 + dof], speed_limit,
                       f"{name} speed DOF {dof}")

    history = result.ptos["relative_heave"]
    _max_error(history.stroke, pto[:, 3], 0.007, "PTO stroke")
    _max_error(history.velocity, pto[:, 9], 0.006, "PTO speed")
    _max_error(history.force, pto[:, 15], 6_500, "PTO total force")
    _max_error(history.force,
               -1_200_000 * history.velocity
               - 50_000 * (history.stroke + 0.1),
               1e-8, "Python PTO force law")
    _max_error(history.absorbed_power,
               1_200_000 * history.velocity**2,
               1e-8, "Python damper dissipation")
    _max_error(pto[:, 15],
               -1_200_000 * pto[:, 9] - 50_000 * (pto[:, 3] + 0.1),
               1e-6, "MATLAB spring and damper force law")
    _max_error(pto[:, 21], pto[:, 15] * pto[:, 9],
               1e-6, "MATLAB spring and damper mechanical power")
    _max_error(-history.force * history.velocity, -pto[:, 21],
               8_000, "PTO mechanical power")

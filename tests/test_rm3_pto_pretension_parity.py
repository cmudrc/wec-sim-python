"""Pair a configured RM3 PTO spring and pretension with pinned WEC-Sim."""

import os
from pathlib import Path

import numpy as np

from wecsim import RegularWave, WEC


ROOT = Path(__file__).resolve().parents[1]
HYDRO = Path(os.environ.get("WEC_SIM_RM3_H5", ROOT / "examples/data/rm3.h5"))
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
PORTABLE = ROOT / "tests/fixtures/matlab/rm3_pto_pretension/source.npz"


def _source():
    """Use fresh MATLAB CSVs when supplied, otherwise the pinned CI trace."""
    if (REFERENCE and (Path(REFERENCE)
                       / "RM3_PTO_PRETENSION_parameters.csv").is_file()):
        source = Path(REFERENCE)
        bodies = [
            np.loadtxt(source / f"RM3_PTO_PRETENSION_Pretension_body{number}.csv",
                       delimiter=",")
            for number in (1, 2)
        ]
        pto = np.loadtxt(
            source / "RM3_PTO_PRETENSION_Pretension_pto1.csv", delimiter=",")
        assert all(values.shape == (4001, 25) for values in (*bodies, pto))
        return {
            "time": pto[:, 0],
            "float_position": bodies[0][:, [1, 3, 5]],
            "float_velocity": bodies[0][:, [7, 9, 11]],
            "spar_position": bodies[1][:, [1, 3, 5]],
            "spar_velocity": bodies[1][:, [7, 9, 11]],
            "pto_stroke": pto[:, 3],
            "pto_velocity": pto[:, 9],
            "pto_force": pto[:, 15],
            "pto_mechanical_power": pto[:, 21],
            "parameters": np.loadtxt(
                source / "RM3_PTO_PRETENSION_parameters.csv", delimiter=",",
                ndmin=1,
            ),
        }
    with np.load(PORTABLE, allow_pickle=False) as saved:
        return {key: saved[key].copy() for key in saved.files}


def _max_error(actual, expected, limit, label):
    actual = np.asarray(actual)
    expected = np.asarray(expected)
    assert actual.shape == expected.shape, label
    assert np.isfinite(actual).all() and np.isfinite(expected).all(), label
    error = float(np.max(np.abs(actual - expected)))
    assert error < limit, f"{label}: {error:.6g} exceeds {limit}"


def test_rm3_configured_spring_and_pretension_against_matlab():
    source = _source()
    np.testing.assert_allclose(source["parameters"], [50_000, 5_000, -0.1],
                               rtol=0, atol=1e-12)
    time = source["time"]
    assert time.shape == (4001,)
    np.testing.assert_allclose(time, np.arange(4001) * .1, rtol=0, atol=1e-10)
    for name in ("float", "spar"):
        assert source[f"{name}_position"].shape == (4001, 3)
        assert source[f"{name}_velocity"].shape == (4001, 3)
    for name in ("pto_stroke", "pto_velocity", "pto_force",
                 "pto_mechanical_power"):
        assert source[name].shape == time.shape
    assert all(np.isfinite(values).all() for values in source.values())
    assert abs(source["pto_force"][0]) > 1_000  # preloaded spring at startup

    wec = WEC("RM3 configured PTO")
    float_body = wec.body("float", HYDRO.resolve(),
                          inertia=(0, 21_306_090.66, 0))
    spar = wec.body("spar", HYDRO.resolve(),
                    inertia=(0, 94_407_091.24, 0))
    wec.floating_joint(
        float_body, spar, damping=1_200_000, stiffness=50_000,
        pretension=5_000,
    )
    result = wec.run(
        RegularWave(2.5, 8), dt=0.1, end_time=400, ramp_time=100,
    )
    assert result.case["pto"]["equilibrium_position"] == -0.1
    _max_error(result.time, time, 1e-10, "time")
    for name, position_limits, speed_limits in (
        ("float", (0.035, 0.006, 0.001), (0.025, 0.005, 0.001)),
        ("spar", (0.045, 0.0025, 0.001), (0.035, 0.002, 0.001)),
    ):
        for column, (dof, position_limit, speed_limit) in enumerate(zip(
            (0, 2, 4), position_limits, speed_limits, strict=True,
        )):
            _max_error(result.bodies[name].position[:, dof],
                       source[f"{name}_position"][:, column], position_limit,
                       f"{name} position DOF {dof}")
            _max_error(result.bodies[name].velocity[:, dof],
                       source[f"{name}_velocity"][:, column], speed_limit,
                       f"{name} speed DOF {dof}")

    history = result.ptos["relative_heave"]
    _max_error(history.stroke, source["pto_stroke"], 0.007, "PTO stroke")
    _max_error(history.velocity, source["pto_velocity"], 0.006, "PTO speed")
    _max_error(history.force, source["pto_force"], 6_500, "PTO total force")
    _max_error(history.force,
               -1_200_000 * history.velocity
               - 50_000 * (history.stroke + 0.1),
               1e-8, "Python PTO force law")
    _max_error(history.absorbed_power,
               1_200_000 * history.velocity**2,
               1e-8, "Python damper dissipation")
    _max_error(source["pto_force"],
               -1_200_000 * source["pto_velocity"]
               - 50_000 * (source["pto_stroke"] + 0.1),
               1e-6, "MATLAB spring and damper force law")
    _max_error(source["pto_mechanical_power"],
               source["pto_force"] * source["pto_velocity"],
               1e-6, "MATLAB spring and damper mechanical power")
    _max_error(-history.force * history.velocity,
               -source["pto_mechanical_power"],
               8_000, "PTO mechanical power")

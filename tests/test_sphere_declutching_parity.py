"""Pair the published Sphere declutching controller with MATLAB dynamics."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim import DeclutchingControl, RegularWave, WEC, WorldPoint


SPHERE_H5 = os.environ.get("WEC_SIM_SPHERE_H5")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
ROOT = Path(__file__).resolve().parents[1]
PORTABLE = Path(__file__).parent / "fixtures/matlab/sphere_declutching.npz"


def _max_error(actual, expected, limit, label):
    assert actual.shape == expected.shape, label
    assert np.isfinite(expected).all(), label
    error = np.max(np.abs(actual - expected))
    assert error < limit, f"{label}: max error {error:.6g} exceeds {limit}"


def test_published_sphere_declutching_against_matlab():
    if SPHERE_H5 and REFERENCE:
        reference = Path(REFERENCE)
        body = np.loadtxt(
            reference / "Sphere_Declutching_Declutching_body1.csv", delimiter=",",
        )
        controller = np.loadtxt(
            reference / "Sphere_Declutching_controller.csv", delimiter=",",
        )
        assert body.shape == (40_001, 25)
        assert controller.shape == (40_001, 13)
        source_time, controller_time = body[:, 0], controller[:, 0]
        source_position, source_velocity = body[:, 3], body[:, 9]
        source_force, source_power = controller[:, 3], controller[:, 9]
        hydro = Path(SPHERE_H5).resolve()
    elif not SPHERE_H5 and not REFERENCE:
        # Ordinary Python CI keeps every published sample used by this gate.
        with np.load(PORTABLE) as source:
            source_time = source["time"]
            controller_time = source_time
            source_position = source["position"]
            source_velocity = source["velocity"]
            source_force = source["force"]
            source_power = source["power"]
        hydro = ROOT / "examples/data/sphere.h5"
    else:
        pytest.fail("provide both WEC_SIM_SPHERE_H5 and WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
    assert source_time.shape == (40_001,)
    wec = WEC("Sphere declutching")
    sphere = wec.body(
        "sphere", hydro,
        inertia=(20_907_301, 21_306_090.66, 37_085_481.11),
    )
    wec.coordinate("heave", sphere.move("heave"))
    wec.pto(
        "declutch", sphere.at(0, 0, 2), WorldPoint(0, 0, 0),
        axis=(0, 0, 1),
        control=DeclutchingControl(gain=232_020, declutch_time=0.8),
    )
    response = wec.run(
        RegularWave(2.5, 9.6664), dt=0.01, end_time=400, ramp_time=100,
    )
    np.testing.assert_allclose(response.time, source_time, rtol=0, atol=1e-10)
    np.testing.assert_allclose(response.time, controller_time, rtol=0, atol=1e-10)

    position = response.bodies["sphere"].position[:, 2]
    velocity = response.bodies["sphere"].velocity[:, 2]
    force = -response.ptos["declutch"].force
    power = -response.ptos["declutch"].absorbed_power
    _max_error(position, source_position, 0.0025, "heave position")
    _max_error(velocity, source_velocity, 0.006, "heave velocity")
    _max_error(power, force * velocity, 1e-8,
               "Python controller force-power relation")
    _max_error(source_power, source_force * source_velocity, 1e-6,
               "MATLAB controller force-power relation")

    # A one-sample shift at a velocity reversal changes the mode and makes
    # a pointwise force maximum misleading. Gate event timing and the force
    # law separately, including power on samples engaged in both runs.
    matlab_on = np.abs(source_force) > 1e-5
    python_on = np.abs(force) > 1e-5
    mismatch_count = np.count_nonzero(matlab_on != python_on)
    assert mismatch_count <= 4, f"{mismatch_count} mismatched controller samples"
    matlab_off = np.flatnonzero(matlab_on[:-1] & ~matlab_on[1:]) + 1
    python_off = np.flatnonzero(python_on[:-1] & ~python_on[1:]) + 1
    assert len(matlab_off) == len(python_off) == 83
    _max_error(response.time[python_off], response.time[matlab_off],
               0.01000001, "declutch event time")
    for on in (matlab_on, python_on):
        off = np.flatnonzero(on[:-1] & ~on[1:]) + 1
        reengage = np.flatnonzero(~on[:-1] & on[1:]) + 1
        _max_error(response.time[reengage[1:]] - response.time[off],
                   np.full(len(off), 0.8), 1e-10,
                   "declutch duration")
    engaged = matlab_on & python_on
    _max_error(force[engaged], source_force[engaged], 1_500,
               "engaged controller force")
    _max_error(power[engaged], source_power[engaged], 1_500,
               "engaged controller power")
    _max_error(force[~matlab_on & ~python_on],
               source_force[~matlab_on & ~python_on], 1e-8,
               "disengaged controller force")
    _max_error(response.raw.pto_generalized_force[:, 0], force, 1e-8,
               "applied generalized controller force")
    python_energy = np.trapezoid(
        response.ptos["declutch"].absorbed_power, response.time,
    )
    matlab_energy = -np.trapezoid(source_power, response.time)
    assert abs(python_energy - matlab_energy) < 15_000, (
        f"absorbed energy differs by {python_energy - matlab_energy:.1f} J"
    )

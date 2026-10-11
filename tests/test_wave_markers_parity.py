"""Pair published RM3 and OSWEC wave-marker runs with MATLAB output."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim import RegularWave, WEC, WorldPoint


REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
pytestmark = pytest.mark.skipif(
    not REFERENCE, reason="paired MATLAB wave-marker output not provided",
)


@pytest.mark.parametrize("family", ["RM3", "OSWEC"])
def test_published_regular_wave_markers(family):
    reference = Path(REFERENCE)
    prefix = f"WAVE_MARKERS_{family}"
    locations = np.loadtxt(reference / f"{prefix}_locations.csv", delimiter=",")
    markers = np.loadtxt(reference / f"{prefix}_markers.csv", delimiter=",")
    origin = np.loadtxt(reference / f"{prefix}_origin.csv", delimiter=",")
    height, period, direction, ramp_time, gravity, depth, source_k = np.loadtxt(
        reference / f"{prefix}_wave_parameters.csv", delimiter=",",
    )

    grid = np.arange(-20, 21, 10)
    x, y = np.meshgrid(grid, grid)
    expected_locations = np.column_stack((x.ravel(order="F"), y.ravel(order="F")))
    np.testing.assert_array_equal(locations, expected_locations)
    assert markers.shape == (4001, 26)
    assert origin.shape == (4001, 2)
    assert np.isfinite(markers).all()
    np.testing.assert_allclose(markers[:, 0], np.arange(4001) * 0.1,
                               rtol=0, atol=1e-10)
    np.testing.assert_allclose(markers[:, 0], origin[:, 0], rtol=0, atol=1e-10)

    wave = RegularWave(height, period, direction)
    actual = wave.elevation_at(
        markers[:, 0], locations, water_depth=depth,
        ramp_time=ramp_time, g=gravity,
    )
    np.testing.assert_allclose(actual, markers[:, 1:], rtol=0, atol=1e-10)
    center = np.flatnonzero(np.all(locations == [0, 0], axis=1))
    assert center.size == 1
    np.testing.assert_allclose(actual[:, center[0]], origin[:, 1],
                               rtol=0, atol=1e-10)
    # Source k is retained to distinguish a dispersion mismatch from a
    # marker-phase mismatch if the pair ever fails.
    omega = 2 * np.pi / period
    residual = omega**2 - gravity * source_k * np.tanh(source_k * depth)
    assert abs(residual) / omega**2 < 1e-10


def _max_error(actual, expected, limit, label):
    assert actual.shape == expected.shape, label
    assert np.isfinite(actual).all() and np.isfinite(expected).all(), label
    error = float(np.max(np.abs(actual - expected)))
    assert error < limit, f"{label}: {error:.6g} exceeds {limit}"


@pytest.mark.parametrize("family", ["RM3", "OSWEC"])
def test_published_wave_marker_body_and_pto_motion(family):
    applications = os.environ.get("WEC_SIM_APPLICATIONS_DIR")
    if not applications:
        pytest.skip("paired MATLAB applications HDF5 not provided")

    hydro = (Path(applications) / "_Common_Input_Files" / family /
             "hydroData" / f"{family.lower()}.h5")
    wec = WEC(f"{family} published Wave Markers")
    if family == "RM3":
        moving = wec.body("float", hydro, inertia=(0, 21_306_090.66, 0))
        other = wec.body("spar", hydro, inertia=(0, 94_407_091.24, 0))
        wec.floating_joint(moving, other, damping=1_200_000,
                           radiation_method="constant")
        body_names = ("float", "spar")
        pto_name = "relative_heave"
        pto_columns = (3, 9, 15, 21)
        pto_limits = (0.006, 0.005, 6000, 7000)
        position_limits = {0: 0.04, 2: 0.006, 4: 0.001}
        speed_limits = {0: 0.03, 2: 0.004, 4: 0.001}
    else:
        moving = wec.body("flap", hydro, mass=127_000,
                          inertia=(0, 1.85e6, 0))
        other = wec.body("base", hydro, mass=999,
                         inertia=(999, 999, 999))
        wec.fixed_hinge(moving, other, location=WorldPoint(0, 0, -10),
                        pto_location=WorldPoint(0, 0, -8.9), damping=12_000)
        body_names = ("flap", "base")
        pto_name = "hinge"
        pto_columns = (5, 11, 17, 23)
        pto_limits = (0.006, 0.005, 60, 30)
        position_limits = {0: 0.03, 2: 0.012, 4: 0.007}
        speed_limits = {0: 0.025, 2: 0.01, 4: 0.006}

    result = wec.run(RegularWave(2.5, 8), dt=0.1, end_time=400,
                     ramp_time=100)
    assert result.time.shape == (4001,)
    reference = Path(REFERENCE)
    prefix = f"WAVE_MARKERS_{family}_{family}"
    for index, name in enumerate(body_names, start=1):
        expected = np.loadtxt(reference / f"{prefix}_body{index}.csv",
                              delimiter=",")
        assert expected.shape == (4001, 25)
        _max_error(result.time, expected[:, 0], 1e-10,
                   f"{family} body {index} time")
        position = result.bodies[name].position
        velocity = result.bodies[name].velocity
        active = index == 1 or family == "RM3"
        for dof in range(6):
            # The published RM3 layout moves in surge, heave, and pitch;
            # the OSWEC base is fixed and both layouts have inactive axes.
            pos_limit = position_limits.get(dof, 1e-10) if active else 1e-10
            vel_limit = speed_limits.get(dof, 1e-10) if active else 1e-10
            _max_error(position[:, dof], expected[:, 1 + dof], pos_limit,
                       f"{family} {name} position DOF {dof}")
            _max_error(velocity[:, dof], expected[:, 7 + dof], vel_limit,
                       f"{family} {name} speed DOF {dof}")

    expected_pto = np.loadtxt(reference / f"{prefix}_pto1.csv", delimiter=",")
    assert expected_pto.shape == (4001, 25)
    pto = result.ptos[pto_name]
    for label, actual, column, limit, sign in (
        ("stroke", pto.stroke, pto_columns[0], pto_limits[0], 1),
        ("speed", pto.velocity, pto_columns[1], pto_limits[1], 1),
        ("force", pto.force, pto_columns[2], pto_limits[2], 1),
        ("absorbed power", pto.absorbed_power, pto_columns[3],
         pto_limits[3], -1),
    ):
        _max_error(actual, sign * expected_pto[:, column], limit,
                   f"{family} PTO {label}")

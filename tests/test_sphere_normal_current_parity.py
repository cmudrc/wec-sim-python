"""Pair a current-loaded normal/tangential Morison Sphere trajectory."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim import Current, RegularWave, WEC
from wecsim.morison import (
    MorisonElement, regular_morison_source_force,
    regular_wave_normal_tangential_morison_terms,
)


SPHERE_H5 = os.environ.get("WEC_SIM_SPHERE_H5")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
pytestmark = pytest.mark.skipif(
    not (SPHERE_H5 and REFERENCE),
    reason="paired normal/tangential current Sphere output and HDF5 absent",
)


def _response(current):
    wec = WEC("Sphere normal/tangential Morison current")
    sphere = wec.body(
        "sphere", Path(SPHERE_H5).resolve(),
        inertia=(20907301, 21306090.66, 37085481.11),
    )
    for axis in ("surge", "heave", "pitch"):
        wec.coordinate(axis, sphere.move(axis))
    wec.morison_element(
        sphere, point=sphere.at(0, 0, -2), axis=(0, 1, 0),
        drag_coefficient=(1, 0, 0),
        added_mass_coefficient=(1, 0, 0),
        area=(100, 0, 0), volume=20,
    )
    return wec.run(
        RegularWave(1, 8, current=current),
        dt=.01, end_time=40, ramp_time=10,
        initial_coordinate={"heave": 1},
    )


def test_normal_current_source_force_and_independent_motion():
    reference = Path(REFERENCE)
    prefix = "SPHERE_MOVING_MORISON_NORMAL"
    no_current = np.loadtxt(
        reference / f"{prefix}_1m-ME_body1.csv", delimiter=",",
    )
    current = np.loadtxt(
        reference / f"{prefix}_CURRENT_1m-ME_body1.csv", delimiter=",",
    )
    direct_delta = np.loadtxt(
        reference / f"{prefix}_CURRENT_direct_delta.csv", delimiter=",",
    )
    wave = np.loadtxt(
        reference / f"{prefix}_CURRENT_wave.csv", delimiter=",",
    )
    assert no_current.shape == current.shape == (4001, 37)
    assert direct_delta.shape == (4001, 7)
    assert wave.shape == (4001, 2)
    assert all(np.isfinite(a).all() for a in (no_current, current, direct_delta, wave))

    element = MorisonElement(
        point=(0, 0, -2), axis=(0, 1, 0),
        drag_coefficient=(1, 0, 0),
        added_mass_coefficient=(1, 0, 0),
        area=(100, 0, 0), volume=20,
    )
    physical_increment = []
    source_increment = []
    for row in current:
        force_args = dict(
            time=row[0], position=row[1:7], velocity=row[7:13],
            wave_height=1, wave_period=8, ramp_time=10,
            water_depth=np.inf, rho=1000,
        )
        current_args = dict(
            current_speed=.8, current_direction=0,
            current_profile="power", current_depth=30,
        )
        physical_increment.append(
            regular_wave_normal_tangential_morison_terms(
                [element], **force_args, **current_args,
            )[0] - regular_wave_normal_tangential_morison_terms(
                [element], **force_args,
            )[0]
        )
        source_args = dict(
            **force_args, acceleration=row[31:37], direction=0,
            body_morison=2, element_axes=[(0, 1, 0)],
        )
        # Replace the physical infinite-depth label with the source's 200 m
        # HDF5 depth. Both use the same deep-water branch at this frequency.
        source_args["water_depth"] = 200
        source_increment.append(
            regular_morison_source_force(
                [element], **source_args, **current_args,
            ) - regular_morison_source_force([element], **source_args)
        )
    physical_increment = np.asarray(physical_increment)
    source_increment = np.asarray(source_increment)
    np.testing.assert_allclose(direct_delta[:, 0], current[:, 0], atol=1e-12)
    assert np.max(np.abs(source_increment - direct_delta[:, 1:])) < 1e-6
    for axis, minimum, label in ((0, 1000, "surge"),
                                 (2, 1000, "heave"),
                                 (4, 10000, "pitch")):
        peak = np.max(np.abs(direct_delta[:, 1 + axis]))
        assert peak > minimum, f"{label} current loading was not exercised"
        error = np.max(np.abs(
            physical_increment[:, axis] - direct_delta[:, 1 + axis]
        ))
        assert error < .02 * peak, (
            f"{label} physical/source current increment differs by {error:.6g}"
        )

    python_current = _response(Current(.8, 0, "power", 30))
    python_no_current = _response(None)
    assert np.max(np.abs(python_current.time - current[:, 0])) < 1e-10
    assert np.max(np.abs(python_current.wave_elevation - wave[:, 1])) < 1e-10
    print("normal/tangential current motion diagnostic:")
    for axis, label in ((0, "surge"), (2, "heave"), (4, "pitch")):
        actual = python_current.bodies["sphere"]
        baseline = python_no_current.bodies["sphere"]
        source_position_increment = (
            current[:, 1 + axis] - no_current[:, 1 + axis]
        )
        python_position_increment = (
            actual.position[:, axis] - baseline.position[:, axis]
        )
        print(label, "absolute position", np.max(np.abs(
            actual.position[:, axis] - current[:, 1 + axis])),
              "absolute speed", np.max(np.abs(
                  actual.velocity[:, axis] - current[:, 7 + axis])),
              "current response peak", np.max(np.abs(source_position_increment)),
              "current response error", np.max(np.abs(
                  python_position_increment - source_position_increment)))

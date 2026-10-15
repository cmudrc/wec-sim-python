"""Pair a sampled-elevation Sphere heave WEC with pinned MATLAB WEC-Sim."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim import ImportedElevationWave, WEC


APPLICATIONS = os.environ.get("WEC_SIM_APPLICATIONS_DIR")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
ROOT = Path(__file__).resolve().parents[1]
PORTABLE = Path(__file__).parent / "fixtures/matlab/sphere_elevation_import"


def _max_error(actual, expected, limit, label):
    actual = np.asarray(actual)
    expected = np.asarray(expected)
    assert actual.shape == expected.shape, label
    assert np.isfinite(actual).all() and np.isfinite(expected).all(), label
    error = float(np.max(np.abs(actual - expected)))
    assert error < limit, f"{label}: max error {error:.6g} exceeds {limit}"


def _source():
    if APPLICATIONS and REFERENCE:
        root = Path(APPLICATIONS)
        reference = Path(REFERENCE)
        body = np.loadtxt(
            reference / "SPHERE_ELEVATION_IMPORT_0m_body1.csv", delimiter=",",
        )
        wave = np.loadtxt(
            reference / "SPHERE_ELEVATION_IMPORT_wave.csv", delimiter=",",
        )
        assert body.shape == (4001, 25)
        assert wave.shape == (4001, 2)
        _max_error(body[:, 0], wave[:, 0], 1e-10, "source time")
        return (body[:, 0], wave[:, 1], body[:, 19:25], body[:, 3],
                body[:, 9], root / "_Common_Input_Files/Sphere/hydroData/sphere.h5",
                root / "Free_Decay/0m/etaData.mat", root)
    if not APPLICATIONS and not REFERENCE:
        with np.load(PORTABLE / "source.npz") as saved:
            channels = tuple(saved[key] for key in
                             ("time", "wave_elevation", "excitation",
                              "position", "velocity"))
        return (*channels, ROOT / "examples/data/sphere.h5",
                PORTABLE / "etaData.mat", ROOT)
    pytest.fail("provide both WEC_SIM_APPLICATIONS_DIR and "
                "WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")


def test_public_sphere_imported_elevation_motion_and_force():
    time, wave, excitation, position, velocity, hydro, eta, root = _source()
    assert time.shape == (4001,)
    assert wave.shape == (4001,)
    assert excitation.shape == (4001, 6)
    wec = WEC("Sphere with imported elevation")
    sphere = wec.body(
        "sphere", hydro,
        inertia=(20_907_301, 21_306_090.66, 37_085_481.11),
    )
    wec.coordinate("heave", sphere.move("heave"))
    result = wec.run(
        ImportedElevationWave(
            eta, reapply_force_ramp=True,
        ),
        dt=0.01, end_time=40, ramp_time=10, radiation_memory=15,
        base_dir=root,
    )
    _max_error(result.time, time, 1e-10, "time")
    _max_error(result.wave_elevation, wave, 1e-12,
               "wave elevation")
    force = dict(result.raw.extra_outputs)["body1_excitation_force"]
    _max_error(force, excitation, 1e-6, "six excitation forces")
    _max_error(result.bodies["sphere"].position[:, 2], position, 1e-4,
               "heave position")
    _max_error(result.bodies["sphere"].velocity[:, 2], velocity, 1e-4,
               "heave velocity")
    assert result.raw.auxiliary_files == (eta.resolve(),)

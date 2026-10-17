"""Pair two-column imported seas with a seeded MATLAB waveClass run."""

import os
from hashlib import sha256
from pathlib import Path

import numpy as np
import pytest
from scipy.io import savemat

from wecsim import ImportedSpectrumWave, WEC
from wecsim.environment.irregularWave import imported_spectrum_components


ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests/fixtures/matlab/imported_spectrum_two_column"
SOURCE = Path(os.environ.get("WEC_SIM_MATLAB_REFERENCE_DIR", FIXTURE))
SPHERE = ROOT / "examples/data/sphere.h5"
RM3 = ROOT / "examples/data/rm3.h5"


def test_pinned_two_column_source_files_are_unchanged():
    expected = {
        "spectrum.mat": "1d79fbb0a4855fddbf9906b7c9e16307a6359c191bee9b512da6964aecba9ea2",
        "bins_phase.csv": "b556a831b75c153739c991b3bd9ef272887bc96a314acde19fda42f9d3b97533",
        "elevation.csv": "328147d68a2d8ad5dcde5afff019cf9170780608f833484187840e03ed3ac5b9",
    }
    for name, digest in expected.items():
        assert sha256((FIXTURE / name).read_bytes()).hexdigest() == digest


def _source(label):
    name = (f"imported_two_column_{label}.csv"
            if SOURCE != FIXTURE else f"{label}.csv")
    return np.loadtxt(SOURCE / name, delimiter=",")


def test_seeded_two_column_import_matches_matlab_bins_phases_and_public_sea():
    """Filter before drawing phases, then carry the sea through WEC.run."""
    bins = _source("bins_phase")
    elevation = _source("elevation")
    components = imported_spectrum_components(
        SPHERE, FIXTURE / "spectrum.mat", seed=1, phase_generator="matlab",
    )
    assert bins.shape == (5, 4) and elevation.shape == (201, 2)
    np.testing.assert_allclose(components.omega, bins[:, 0], rtol=0, atol=1e-13)
    np.testing.assert_allclose(components.d_omega, bins[:, 1], rtol=0, atol=1e-13)
    np.testing.assert_allclose(components.spectral_amplitude / 2, bins[:, 2],
                               rtol=0, atol=1e-13)
    np.testing.assert_allclose(components.phase[:, 0], bins[:, 3],
                               rtol=0, atol=1e-13)

    wec = WEC("seeded imported Sphere sea")
    sphere = wec.body("sphere", SPHERE)
    wec.coordinate("heave", sphere.move("heave"))
    response = wec.run(
        ImportedSpectrumWave(FIXTURE / "spectrum.mat", seed=1,
                             phase_generator="matlab"),
        dt=.1, end_time=20, ramp_time=1, radiation_memory=15,
    )
    assert response.case["wave"]["seed"] == 1
    assert response.case["wave"]["phase_generator"] == "matlab"
    np.testing.assert_allclose(response.time, elevation[:, 0], rtol=0,
                               atol=1e-13)
    np.testing.assert_allclose(response.wave_elevation, elevation[:, 1],
                               rtol=0, atol=1e-12)
    assert np.ptp(response.wave_elevation) > .5
    assert np.isfinite(response.bodies["sphere"].position).all()


def test_two_column_import_reaches_floating_joint_with_the_same_sea():
    """The RM3 joint and mapped-coordinate adapters share phase semantics."""
    wec = WEC("seeded imported RM3 sea")
    float_body = wec.body("float", RM3, inertia=(0, 21_306_090.66, 0))
    spar = wec.body("spar", RM3, inertia=(0, 94_407_091.24, 0))
    wec.floating_joint(float_body, spar, damping=1_200_000)
    response = wec.run(
        ImportedSpectrumWave(FIXTURE / "spectrum.mat", seed=1,
                             phase_generator="matlab"),
        dt=.1, end_time=2, ramp_time=1, radiation_memory=1,
    )
    elevation = _source("elevation")
    np.testing.assert_allclose(response.wave_elevation, elevation[:21, 1],
                               rtol=0, atol=1e-12)
    assert np.isfinite(response.bodies["float"].position).all()


def test_imported_spectrum_phase_selection_is_explicit(tmp_path):
    values = np.array([[.1, .2], [.15, .4], [.2, .2]])
    path = tmp_path / "spectrum.mat"
    savemat(path, {"spectrumData": values})
    first = imported_spectrum_components(SPHERE, path, seed=1)
    again = imported_spectrum_components(SPHERE, path, seed=1)
    second = imported_spectrum_components(SPHERE, path, seed=2)
    np.testing.assert_array_equal(first.phase, again.phase)
    assert not np.array_equal(first.phase, second.phase)
    with pytest.raises(ValueError, match="positive integer substream seed"):
        imported_spectrum_components(SPHERE, path, phase_generator="matlab")

    savemat(path, {"spectrumData": np.column_stack((values, [0, 1, 2]))})
    with pytest.raises(ValueError, match="cannot be overridden"):
        imported_spectrum_components(SPHERE, path, seed=1)
    with pytest.raises(ValueError, match="phase_generator"):
        ImportedSpectrumWave(path, phase_generator="unknown").as_case()

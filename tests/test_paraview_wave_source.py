"""Pair Python VTP wave meshes with the pinned MATLAB ParaView writer."""

import os
from pathlib import Path
from xml.etree import ElementTree

import h5py
import numpy as np
import pytest

from wecsim.visualization.paraviewClass import ParaviewClass
from wecsim.environment.waveClass import WaveClass


REFERENCE = os.environ.get("WEC_SIM_MATLAB_PARAVIEW_WAVE_DIR")
pytestmark = pytest.mark.skipif(
    not REFERENCE, reason="pinned MATLAB ParaView wave grids not provided",
)


def _mesh(path):
    piece = ElementTree.parse(path).find("./PolyData/Piece")
    points = np.fromstring(
        piece.findtext("./Points/DataArray"), sep=" ",
    ).reshape(-1, 3)
    connectivity = np.fromstring(
        piece.findtext("./Polys/DataArray[@Name='connectivity']"),
        sep=" ", dtype=int,
    ).reshape(-1, 4)
    offsets = np.fromstring(
        piece.findtext("./Polys/DataArray[@Name='offsets']"),
        sep=" ", dtype=int,
    )
    assert int(piece.get("NumberOfPoints")) == len(points)
    assert int(piece.get("NumberOfPolys")) == len(connectivity)
    return points, connectivity, offsets


def test_regular_wave_vtp_against_pinned_matlab(tmp_path):
    reference = Path(REFERENCE)
    amplitude, wavenumber, omega, direction, depth = np.loadtxt(
        reference / "parameters.csv", delimiter=",",
    )
    wave = WaveClass("regular")
    wave.A = amplitude
    wave.k = wavenumber
    wave.w = omega
    wave.waveDir = [direction]
    wave.waterDepth = depth
    paths = ParaviewClass(wave).write_paraview_vtp_wave(
        [2, 3], tmp_path, domain_size=6, num_points_x=3, num_points_y=2,
    )
    assert (tmp_path / "ground.txt").read_text() == (
        reference / "ground.txt"
    ).read_text()
    for index, path in enumerate(paths, start=1):
        source = reference / "waves" / f"waves_{index}.vtp"
        actual_points, actual_cells, actual_offsets = _mesh(path)
        source_points, source_cells, source_offsets = _mesh(source)
        np.testing.assert_allclose(actual_points, source_points, rtol=0, atol=1e-5)
        np.testing.assert_array_equal(actual_cells, source_cells)
        np.testing.assert_array_equal(actual_offsets, source_offsets)


def test_imported_elevation_vtp_against_pinned_matlab(tmp_path):
    reference = Path(REFERENCE) / "imported"
    record = np.loadtxt(reference / "record.csv", delimiter=",")
    assert record.shape == (16, 2)
    sample_time = np.arange(16) * .25
    expected_elevation = (
        .7 * np.cos(np.pi * sample_time + .35)
        + .2 * np.cos(2 * np.pi * sample_time - .8)
    )
    np.testing.assert_allclose(record[:, 0], sample_time, rtol=0, atol=1e-12)
    np.testing.assert_allclose(record[:, 1], expected_elevation, rtol=0, atol=1e-12)

    wave = WaveClass("etaImport")
    wave.waveAmpTime = [record[:, 0], record[:, 1]]
    wave.waveDir = [60]
    wave.waterDepth = 30
    path, = ParaviewClass(wave).write_paraview_vtp_wave(
        [.375], tmp_path, domain_size=2,
        num_points_x=3, num_points_y=2, g=9.7,
    )
    assert (tmp_path / "ground.txt").read_text() == (
        reference / "ground.txt"
    ).read_text()
    actual_points, actual_cells, actual_offsets = _mesh(path)
    source_points, source_cells, source_offsets = _mesh(
        reference / "waves/waves_1.vtp"
    )
    np.testing.assert_allclose(actual_points, source_points, rtol=0, atol=1e-5)
    np.testing.assert_array_equal(actual_cells, source_cells)
    np.testing.assert_array_equal(actual_offsets, source_offsets)


def test_published_rm3_imported_wave_vtp_against_pinned_matlab(tmp_path):
    """The published MooringMatrix sea also produces the source wave surface."""
    reference = Path(REFERENCE) / "rm3_mooring_matrix"
    wave = WaveClass("etaImport")
    fixture = Path(__file__).parent / "fixtures/matlab/rm3_mooring_matrix"
    wave.etaDataFile = str(fixture / "etaData.mat")
    with h5py.File(fixture / "rm3.h5") as hydro:
        bem_frequency = hydro["simulation_parameters/w"][:].ravel()
        assert hydro["simulation_parameters/water_depth"][()] == b"infinite"
    wave.waveSetup(bem_frequency, "infinite", 40, .01, 40000, 9.81, 1000, 400)
    assert wave.waterDepth == 200

    source_record = np.loadtxt(reference / "wave_amp_time.csv", delimiter=",")
    np.testing.assert_allclose(
        np.asarray(wave.waveAmpTime).T, source_record, rtol=0, atol=1e-10,
    )
    frames = [0, 80, 240, 400]
    paths = ParaviewClass(wave).write_paraview_vtp_wave(
        frames, tmp_path, domain_size=4,
        num_points_x=4, num_points_y=3,
    )
    assert (tmp_path / "ground.txt").read_text() == (
        reference / "ground.txt"
    ).read_text()
    for path, index in zip(paths, range(1, len(frames) + 1), strict=True):
        actual_points, actual_cells, actual_offsets = _mesh(path)
        source_points, source_cells, source_offsets = _mesh(
            reference / f"waves/waves_{index}.vtp"
        )
        np.testing.assert_allclose(actual_points, source_points, rtol=0, atol=1e-5)
        np.testing.assert_array_equal(actual_cells, source_cells)
        np.testing.assert_array_equal(actual_offsets, source_offsets)

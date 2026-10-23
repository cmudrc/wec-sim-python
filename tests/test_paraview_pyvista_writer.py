"""Check the optional VTP writer against the source-compatible serializer."""

import os
from pathlib import Path
from xml.etree import ElementTree

import numpy as np
import pytest

from wecsim.visualization.paraviewClass import ParaviewClass
from wecsim.environment.waveClass import WaveClass


pyvista = pytest.importorskip("pyvista")


def _polydata(path):
    piece, = ElementTree.parse(path).findall("./PolyData/Piece")
    points = np.fromstring(piece.findtext("./Points/DataArray"), sep=" ").reshape(-1, 3)
    polys = piece.find("./Polys")
    cells = np.fromstring(
        polys.find("./DataArray[@Name='connectivity']").text,
        sep=" ", dtype=int,
    )
    offsets = np.fromstring(
        polys.find("./DataArray[@Name='offsets']").text,
        sep=" ", dtype=int,
    )
    fields = {
        array.get("Name"): np.fromstring(array.text, sep=" ")
        for array in piece.findall("./CellData/DataArray")
    }
    return points, cells, offsets, fields


def _compare_frames(source, optional, *, point_limit=1e-5, field_limit=1e-6):
    assert [path.name for path in source] == [path.name for path in optional]
    for source_path, optional_path in zip(source, optional):
        source_points, source_cells, source_offsets, source_fields = _polydata(source_path)
        points, cells, offsets, fields = _polydata(optional_path)
        np.testing.assert_allclose(points, source_points, rtol=0, atol=point_limit)
        np.testing.assert_array_equal(cells, source_cells)
        np.testing.assert_array_equal(offsets, source_offsets)
        assert fields.keys() == source_fields.keys()
        for name in fields:
            np.testing.assert_allclose(
                fields[name], source_fields[name], rtol=0, atol=field_limit,
            )
        mesh = pyvista.read(optional_path)
        assert mesh.n_points == len(points)
        assert mesh.n_cells == len(offsets)
        assert set(mesh.cell_data.keys()) == set(fields)


def test_wave_vtp_pyvista_matches_source_layout_and_values(tmp_path):
    class RegularWave:
        wType = "regular"
        waveDir = [30]
        A = 1.25
        k = 0.3
        w = 0.8
        waterDepth = 30

    writer = ParaviewClass(RegularWave())
    kwargs = dict(domain_size=6, num_points_x=4, num_points_y=3)
    source_dir = tmp_path / "source"
    optional_dir = tmp_path / "pyvista"
    source = writer.write_paraview_vtp_wave([2, 3], source_dir, **kwargs)
    optional = writer.write_paraview_vtp_wave(
        [2, 3], optional_dir, writer="pyvista", **kwargs,
    )
    assert (source_dir / "ground.txt").read_bytes() == (
        optional_dir / "ground.txt"
    ).read_bytes()
    _compare_frames(source, optional)


def test_body_vtp_pyvista_matches_source_layout_and_pressure(tmp_path):
    writer = ParaviewClass(None)
    kwargs = dict(
        body_name="flap",
        vertices=np.array([[0, 0, 0], [2, 0, 0], [0, 1, 0], [0, 0, 2]]),
        faces=np.array([[0, 1, 2], [0, 2, 3]]),
        poses=np.array([[0, 0, 0, 0, 0, 0], [1, -2, 3, .2, -.3, .4]]),
        cell_areas=[1, 1],
        hydrostatic_pressure=[[10, 20], [30, 40]],
        nonlinear_wave_pressure=[[1, 2], [3, 4]],
        linear_wave_pressure=[[5, 6], [7, 8]],
    )
    source = writer.write_paraview_vtp([2, 3], tmp_path / "source", **kwargs)
    optional = writer.write_paraview_vtp(
        [2, 3], tmp_path / "pyvista", writer="pyvista", **kwargs,
    )
    _compare_frames(source, optional)


def test_invalid_writer_rejected(tmp_path):
    class FlatWave:
        wType = "noWave"
        waterDepth = 10

    with pytest.raises(ValueError, match="writer"):
        ParaviewClass(FlatWave()).write_paraview_vtp_wave(
            [0], tmp_path, domain_size=1, num_points_x=2, num_points_y=2,
            writer="unknown",
        )


@pytest.mark.skipif(
    not os.environ.get("WEC_SIM_MATLAB_PARAVIEW_WAVE_DIR"),
    reason="pinned MATLAB wave VTP output not provided",
)
def test_pyvista_wave_against_pinned_matlab(tmp_path):
    reference = Path(os.environ["WEC_SIM_MATLAB_PARAVIEW_WAVE_DIR"])
    amplitude, wavenumber, omega, direction, depth = np.loadtxt(
        reference / "parameters.csv", delimiter=",",
    )
    wave = WaveClass("regular")
    wave.A, wave.k, wave.w = amplitude, wavenumber, omega
    wave.waveDir = [direction]
    wave.waterDepth = depth
    paths = ParaviewClass(wave).write_paraview_vtp_wave(
        [2, 3], tmp_path, domain_size=6, num_points_x=3,
        num_points_y=2, writer="pyvista",
    )
    _compare_frames(
        [reference / "waves" / f"waves_{index}.vtp" for index in (1, 2)],
        paths,
    )


@pytest.mark.skipif(
    not os.environ.get("WEC_SIM_MATLAB_PARAVIEW_BODY_DIR"),
    reason="pinned MATLAB body VTP output not provided",
)
def test_pyvista_body_against_pinned_matlab(tmp_path):
    reference = Path(os.environ["WEC_SIM_MATLAB_PARAVIEW_BODY_DIR"])
    paths = ParaviewClass(None).write_paraview_vtp(
        [2, 3], tmp_path, body_name="flap",
        vertices=[[0, 0, 0], [2, 0, 0], [0, 1, 0], [0, 0, 2]],
        faces=[[0, 1, 2], [0, 2, 3]],
        poses=[[0, 0, 0, 0, 0, 0], [1, -2, 3, .2, -.3, .4]],
        cell_areas=[1, 1],
        hydrostatic_pressure=[[10, 20], [30, 40]],
        nonlinear_wave_pressure=[[1, 2], [3, 4]],
        linear_wave_pressure=[[5, 6], [7, 8]],
        writer="pyvista",
    )
    _compare_frames(
        [reference / "body1_flap" / f"flap_{index}.vtp" for index in (1, 2)],
        paths,
    )

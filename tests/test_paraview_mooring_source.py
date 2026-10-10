"""Pair MoorDyn line VTP output with the pinned MATLAB writer."""

import os
from pathlib import Path
from xml.etree import ElementTree

import numpy as np
import pytest

from wecsim.paraviewClass import ParaviewClass


REFERENCE = os.environ.get("WEC_SIM_MATLAB_PARAVIEW_MOORING_DIR")
SOURCE_TIMES = np.array([0., 1., 2.])
FRAME_TIMES = np.array([.25, 1.25])


def _inputs():
    lines = []
    tensions = []
    for line, node_count in ((1, 3), (2, 4)):
        nodes = np.empty((len(SOURCE_TIMES), node_count, 3))
        for node in range(node_count):
            nodes[:, node, 0] = 10 * line + node + .4 * SOURCE_TIMES
            nodes[:, node, 1] = 2 * line - .25 * node + .3 * SOURCE_TIMES
            nodes[:, node, 2] = -line - 1.5 * node + .2 * SOURCE_TIMES
        force = np.empty((len(SOURCE_TIMES), node_count - 1))
        for segment in range(1, node_count):
            force[:, segment - 1] = (
                100 * line + 10 * segment + np.array([0, 2.5, 7])
            )
        lines.append(nodes)
        tensions.append(force)
    return lines, tensions


def _read(path):
    pieces = ElementTree.parse(path).findall("./PolyData/Piece")
    result = []
    for piece in pieces:
        points = np.fromstring(
            piece.findtext("./Points/DataArray"), sep=" ",
        ).reshape(-1, 3)
        connectivity = np.fromstring(
            piece.findtext("./Lines/DataArray[@Name='connectivity']"),
            sep=" ", dtype=int,
        ).reshape(-1, 2)
        offsets = np.fromstring(
            piece.findtext("./Lines/DataArray[@Name='offsets']"),
            sep=" ", dtype=int,
        )
        tension = np.fromstring(
            piece.findtext("./CellData/DataArray[@Name='Segment Tension']"),
            sep=" ",
        )
        assert int(piece.get("NumberOfPoints")) == len(points)
        assert int(piece.get("NumberOfLines")) == len(connectivity)
        result.append((points, connectivity, offsets, tension))
    return result


def test_mooring_vtp_interpolates_nodes_and_tensions(tmp_path):
    lines, tensions = _inputs()
    paths = ParaviewClass(None).write_paraview_vtp_mooring(
        SOURCE_TIMES, FRAME_TIMES, tmp_path,
        lines=lines, tensions=tensions,
    )
    assert [path.name for path in paths] == ["mooring_1.vtp", "mooring_2.vtp"]
    for frame, path in enumerate(paths):
        pieces = _read(path)
        assert len(pieces) == 2
        for line, (points, cells, offsets, load) in enumerate(pieces):
            expected = np.array([
                [np.interp(FRAME_TIMES[frame], SOURCE_TIMES,
                           lines[line][:, node, axis])
                 for axis in range(3)]
                for node in range(lines[line].shape[1])
            ])
            np.testing.assert_allclose(points, expected, rtol=0, atol=1e-5)
            np.testing.assert_array_equal(
                cells, np.column_stack((
                    np.arange(len(points) - 1), np.arange(1, len(points)),
                )),
            )
            np.testing.assert_array_equal(
                offsets, 2 * np.arange(1, len(points)),
            )
            np.testing.assert_allclose(
                load,
                [np.interp(FRAME_TIMES[frame], SOURCE_TIMES,
                           tensions[line][:, segment])
                 for segment in range(len(points) - 1)],
                rtol=0, atol=1e-6,
            )
    with pytest.raises(ValueError, match="in-range times"):
        ParaviewClass(None).write_paraview_vtp_mooring(
            SOURCE_TIMES, [-.1], tmp_path, lines=lines, tensions=tensions,
        )


@pytest.mark.skipif(not REFERENCE,
                    reason="pinned MATLAB mooring VTP output not provided")
def test_mooring_vtp_against_pinned_matlab(tmp_path):
    lines, tensions = _inputs()
    paths = ParaviewClass(None).write_paraview_vtp_mooring(
        SOURCE_TIMES, FRAME_TIMES, tmp_path,
        lines=lines, tensions=tensions,
    )
    for frame, actual_path in enumerate(paths, start=1):
        expected_path = (Path(REFERENCE) / "mooring1"
                         / f"mooring_{frame}.vtp")
        actual = _read(actual_path)
        expected = _read(expected_path)
        assert len(actual) == len(expected) == 2
        for (points, cells, offsets, force), source in zip(actual, expected):
            source_points, source_cells, source_offsets, source_force = source
            np.testing.assert_allclose(points, source_points,
                                       rtol=0, atol=1e-5)
            np.testing.assert_array_equal(cells, source_cells)
            np.testing.assert_array_equal(offsets, source_offsets)
            np.testing.assert_allclose(force, source_force,
                                       rtol=0, atol=1e-6)

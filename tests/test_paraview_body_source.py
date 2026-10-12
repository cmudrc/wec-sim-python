"""Check body VTP meshes against the pinned MATLAB ParaView writer."""

import os
from pathlib import Path
from xml.etree import ElementTree

import numpy as np
import pytest

from wecsim.visualization.paraviewClass import ParaviewClass


REFERENCE = os.environ.get("WEC_SIM_MATLAB_PARAVIEW_BODY_DIR")
VERTICES = np.array([[0, 0, 0], [2, 0, 0], [0, 1, 0], [0, 0, 2]])
FACES = np.array([[0, 1, 2], [0, 2, 3]])
POSES = np.array([[0, 0, 0, 0, 0, 0], [1, -2, 3, .2, -.3, .4]])
PRESSURES = {
    "hydrostatic_pressure": np.array([[10, 20], [30, 40]]),
    "nonlinear_wave_pressure": np.array([[1, 2], [3, 4]]),
    "linear_wave_pressure": np.array([[5, 6], [7, 8]]),
}


def _mesh(path):
    piece = ElementTree.parse(path).find("./PolyData/Piece")
    points = np.fromstring(
        piece.findtext("./Points/DataArray"), sep=" ",
    ).reshape(-1, 3)
    faces = np.fromstring(
        piece.findtext("./Polys/DataArray[@Name='connectivity']"),
        sep=" ", dtype=int,
    ).reshape(-1, 3)
    offsets = np.fromstring(
        piece.findtext("./Polys/DataArray[@Name='offsets']"),
        sep=" ", dtype=int,
    )
    fields = {
        array.get("Name"): np.fromstring(array.text, sep=" ")
        for array in piece.findall("./CellData/DataArray")
    }
    assert int(piece.get("NumberOfPoints")) == len(points)
    assert int(piece.get("NumberOfPolys")) == len(faces)
    return points, faces, offsets, fields


def _write(directory):
    return ParaviewClass(None).write_paraview_vtp(
        [2, 3], directory, body_name="flap", vertices=VERTICES,
        faces=FACES, poses=POSES, cell_areas=[1, 1], **PRESSURES,
    )


def test_body_vtp_has_mesh_and_pressure_fields(tmp_path):
    paths = _write(tmp_path)
    assert [path.name for path in paths] == ["flap_1.vtp", "flap_2.vtp"]
    first, faces, offsets, fields = _mesh(paths[0])
    np.testing.assert_array_equal(first, VERTICES)
    np.testing.assert_array_equal(faces, FACES)
    np.testing.assert_array_equal(offsets, [3, 6])
    np.testing.assert_array_equal(fields["CellArea"], [1, 1])
    np.testing.assert_array_equal(fields["HydrostaticPressure"], [10, 20])
    np.testing.assert_array_equal(fields["WavePressureNonLinear"], [1, 2])
    np.testing.assert_array_equal(fields["WavePressureLinear"], [5, 6])
    with pytest.raises(ValueError, match="triangles"):
        ParaviewClass(None).write_paraview_vtp(
            [2], tmp_path, body_name="flap", vertices=VERTICES,
            faces=np.array([[0, 1, 4]]), poses=POSES[:1],
        )


@pytest.mark.skipif(not REFERENCE,
                    reason="pinned MATLAB body VTP output not provided")
def test_body_vtp_against_pinned_matlab(tmp_path):
    reference = Path(REFERENCE)
    for index, path in enumerate(_write(tmp_path), start=1):
        source = reference / "body1_flap" / f"flap_{index}.vtp"
        actual_points, actual_faces, actual_offsets, actual_fields = _mesh(path)
        source_points, source_faces, source_offsets, source_fields = _mesh(source)
        np.testing.assert_allclose(actual_points, source_points, rtol=0, atol=1e-5)
        np.testing.assert_array_equal(actual_faces, source_faces)
        np.testing.assert_array_equal(actual_offsets, source_offsets)
        assert actual_fields.keys() == source_fields.keys()
        for name in actual_fields:
            np.testing.assert_allclose(actual_fields[name], source_fields[name],
                                       rtol=0, atol=1e-6)

"""Pair every published OSWEC nonlinear ParaView wave and body frame."""

import os
from pathlib import Path
from xml.etree import ElementTree

import numpy as np
import pytest
from scipy.io import loadmat
from scipy.spatial import cKDTree
import trimesh

from wecsim.paraviewClass import ParaviewClass
from wecsim.waveClass import WaveClass


REFERENCE = os.environ.get("WEC_SIM_MATLAB_OSWEC_PRESSURE_DIR")
VTP = os.environ.get("WEC_SIM_OSWEC_VTP_DIR")
FLAP = os.environ.get("WEC_SIM_OSWEC_FLAP_STL")
BASE = os.environ.get("WEC_SIM_OSWEC_BASE_STL")
pytestmark = pytest.mark.skipif(
    not (REFERENCE and VTP and FLAP and BASE),
    reason="pinned published OSWEC VTP, poses, and meshes not provided",
)


def _mesh(path):
    piece = ElementTree.parse(path).find("./PolyData/Piece")
    points = np.fromstring(piece.findtext("./Points/DataArray"), sep=" ").reshape(-1, 3)
    connectivity = np.fromstring(
        piece.findtext("./Polys/DataArray[@Name='connectivity']"),
        sep=" ", dtype=int,
    )
    offsets = np.fromstring(
        piece.findtext("./Polys/DataArray[@Name='offsets']"),
        sep=" ", dtype=int,
    )
    fields = {
        item.get("Name"): np.fromstring(item.text, sep=" ")
        for item in piece.findall("./CellData/DataArray")
    }
    assert len(points) == int(piece.get("NumberOfPoints"))
    return points, connectivity, offsets, fields


def _source():
    return loadmat(Path(REFERENCE) / "source.mat")


def _time(source):
    times = source["bodyTime"].ravel()
    np.testing.assert_allclose(times, np.arange(1201) * .1, rtol=0, atol=1e-8)
    return times


def test_published_oswec_wave_vtp_history(tmp_path):
    """Check every wave point, quad, and ground parameter on source times."""
    source = _source()
    times = _time(source)
    domain = float(source["domainSize"].item())
    depth = float(source["waterDepth"].item())
    nx = int(source["numPointsX"].item())
    ny = int(source["numPointsY"].item())
    moorings = int(source["numMoorDyn"].item())
    assert (domain, depth, nx, ny, moorings) == (200, 10.9, 50, 50, 0)
    wave = WaveClass("regular")
    wave.A = 1.25
    wave.w = 2 * np.pi / 8
    wave.k = float(source["waveNumber"].item())
    wave.waveDir = source["waveDirection"].ravel()
    wave.waterDepth = depth
    paths = ParaviewClass(wave).write_paraview_vtp_wave(
        times, tmp_path, domain_size=domain,
        num_points_x=nx, num_points_y=ny, num_moordyn=moorings,
    )
    assert len(paths) == 1201
    reference = Path(VTP)
    assert (tmp_path / "ground.txt").read_text() == (
        reference / "ground.txt"
    ).read_text()
    max_point_error = 0.0
    for frame, path in enumerate(paths, start=1):
        points, quads, offsets, fields = _mesh(path)
        expected, source_quads, source_offsets, source_fields = _mesh(
            reference / "waves" / f"waves_{frame}.vtp"
        )
        assert points.shape == expected.shape == (2500, 3)
        assert fields == source_fields == {}
        max_point_error = max(
            max_point_error, float(np.max(np.abs(points - expected))),
        )
        np.testing.assert_array_equal(quads, source_quads)
        np.testing.assert_array_equal(offsets, source_offsets)
    print(f"OSWEC wave VTP point maximum difference: {max_point_error:.8g} m")
    assert max_point_error < 1.01e-5


@pytest.mark.parametrize(
    ("body_index", "name", "stl", "pose_key", "time_key"),
    (
        (1, "flap", FLAP, "pose", "bodyTime"),
        (2, "base", BASE, "basePose", "baseTime"),
    ),
)
def test_published_oswec_body_vtp_history(
    tmp_path, body_index, name, stl, pose_key, time_key,
):
    """Check all facets and cell fields using saved source body histories."""
    source = _source()
    times = _time(source)
    np.testing.assert_allclose(source[time_key].ravel(), times, rtol=0, atol=1e-8)
    poses = source[pose_key]
    assert poses.shape == (1201, 6)
    mesh = trimesh.load_mesh(stl, process=False)
    kwargs = {}
    if body_index == 1:
        assert len(mesh.faces) == 1042
        np.testing.assert_allclose(mesh.triangles_center, source["centers"],
                                   rtol=0, atol=1e-9)
        kwargs = dict(
            cell_areas=source["areas"].ravel(),
            hydrostatic_pressure=source["hydrostatic"],
            nonlinear_wave_pressure=source["nonlinearWave"],
            linear_wave_pressure=source["linearWave"],
        )
    paths = ParaviewClass(None).write_paraview_vtp(
        times, tmp_path, body_name=name, body_index=body_index,
        vertices=mesh.vertices, faces=mesh.faces, poses=poses, **kwargs,
    )
    assert len(paths) == 1201
    max_center_error = 0.0
    max_vertex_error = 0.0
    max_area_error = 0.0
    max_field_error = {}
    for frame, path in enumerate(paths, start=1):
        points, connectivity, offsets, fields = _mesh(path)
        expected, source_connectivity, source_offsets, source_fields = _mesh(
            Path(VTP) / f"body{body_index}_{name}" / f"{name}_{frame}.vtp"
        )
        assert points.shape == expected.shape
        assert fields.keys() == source_fields.keys()
        assert set(fields) == (
            {"CellArea", "HydrostaticPressure", "WavePressureNonLinear",
             "WavePressureLinear"} if body_index == 1 else {"CellArea"}
        )
        faces = connectivity.reshape(-1, 3)
        source_faces = source_connectivity.reshape(-1, 3)
        np.testing.assert_array_equal(offsets, source_offsets)
        actual_triangles = points[faces]
        source_triangles = expected[source_faces]
        assert actual_triangles.shape == source_triangles.shape
        actual_center = actual_triangles.mean(axis=1)
        source_center = source_triangles.mean(axis=1)
        distance, match = cKDTree(actual_center).query(source_center)
        assert len(np.unique(match)) == len(match)
        max_center_error = max(max_center_error, float(np.max(distance)))
        vertex_distance = np.linalg.norm(
            source_triangles[:, :, None, :]
            - actual_triangles[match][:, None, :, :], axis=3,
        )
        max_vertex_error = max(
            max_vertex_error,
            float(np.max(np.min(vertex_distance, axis=2))),
            float(np.max(np.min(vertex_distance, axis=1))),
        )
        actual_area = np.linalg.norm(np.cross(
            actual_triangles[:, 1] - actual_triangles[:, 0],
            actual_triangles[:, 2] - actual_triangles[:, 0],
        ), axis=1) / 2
        source_area = np.linalg.norm(np.cross(
            source_triangles[:, 1] - source_triangles[:, 0],
            source_triangles[:, 2] - source_triangles[:, 0],
        ), axis=1) / 2
        max_area_error = max(max_area_error, float(np.max(
            np.abs(actual_area[match] - source_area),
        )))
        for field, values in fields.items():
            expected_values = source_fields[field]
            assert len(values) == len(expected_values) == len(faces)
            error = float(np.max(np.abs(values[match] - expected_values)))
            max_field_error[field] = max(max_field_error.get(field, 0), error)
            limit = 2e-6 * np.abs(expected_values) + 1e-4
            assert np.all(np.abs(values[match] - expected_values) < limit), (
                name, frame, field, error,
            )
    print(f"OSWEC {name} VTP facet-center difference: {max_center_error:.8g} m")
    print(f"OSWEC {name} VTP facet-vertex difference: {max_vertex_error:.8g} m")
    print(f"OSWEC {name} VTP facet-area difference: {max_area_error:.8g} m²")
    print(f"OSWEC {name} VTP maximum fields: {max_field_error}")
    assert max_center_error < 3e-5
    assert max_vertex_error < 3e-5
    assert max_area_error < 1e-3

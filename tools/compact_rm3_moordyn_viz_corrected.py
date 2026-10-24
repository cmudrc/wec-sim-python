"""Compact the pinned RM3 MoorDyn corrected source and its actual VTP scene.

Input directories are artifacts of the MATLAB generators in ``tests/matlab``.
The original 0.01 s fairlead grid and all published 0.1 s scene frames are
retained in the portable fixture.
"""

from argparse import ArgumentParser
from pathlib import Path
from xml.etree import ElementTree

import numpy as np


KINDS = {
    "phase": (1000, 1),
    "wave": (801, 2),
    "body1": (801, 25),
    "body2": (801, 25),
    "pto": (801, 25),
    "mooring": (801, 19),
    "fairlead": (801, 4),
    "line1": (801, 84),
    "line2": (801, 84),
    "line3": (801, 84),
}
WAVE_SAMPLE_INDICES = np.array([0, 499, 999, 1000, 1500, 1999])


def _read(path: Path, shape: tuple[int, int]) -> np.ndarray:
    values = np.loadtxt(path, delimiter=",", ndmin=2)
    if values.shape != shape or not np.isfinite(values).all():
        raise ValueError(f"unexpected source data in {path}: {values.shape}")
    return values


def _source(root: Path, prefix: str = "", *, stride: int = 1) -> dict[str, np.ndarray]:
    records = {}
    for kind, shape in KINDS.items():
        full_shape = (8001, shape[1]) if stride == 10 and kind not in (
            "phase", "line1", "line2", "line3"
        ) else shape
        data = _read(root / f"{prefix}{kind}.csv", full_shape)
        records[kind] = data[::stride] if full_shape[0] == 8001 else data
        if kind == "fairlead" and stride == 10:
            records["fairlead_dense"] = data
    if stride == 1:
        records["fairlead_dense"] = _read(root / "fairlead_dense.csv", (8001, 4))
    np.testing.assert_allclose(records["fairlead_dense"][::10],
                               records["fairlead"], rtol=0, atol=1e-8)
    return records


def _points(root: ElementTree.Element, selector: str) -> np.ndarray:
    data = root.findtext(selector)
    if data is None:
        raise ValueError(f"VTP missing {selector}")
    points = np.fromstring(data, sep=" ")
    if points.size % 3 or not np.isfinite(points).all():
        raise ValueError(f"invalid VTP points in {selector}")
    return points.reshape(-1, 3)


def _scene(scene: Path) -> dict[str, np.ndarray]:
    wave_samples = np.empty((801, len(WAVE_SAMPLE_INDICES), 3))
    wave_extrema = np.empty((801, 2))
    body_stats = np.empty((2, 801, 11))
    line_nodes = np.empty((3, 801, 21, 3))
    line_tension = np.empty((3, 801, 20))
    if (scene / "ground.txt").read_text() != "300\n70\n1\n":
        raise ValueError("corrected scene ground settings changed")
    for frame in range(801):
        number = frame + 1
        wave = ElementTree.parse(
            scene / "waves" / f"waves_{number}.vtp"
        ).getroot()
        points = _points(wave, "./PolyData/Piece/Points/DataArray")
        wave_piece = wave.find("./PolyData/Piece")
        wave_cells = np.fromstring(wave_piece.findtext(
            "./Polys/DataArray[@Name='connectivity']"
        ), sep=" ", dtype=int)
        wave_offsets = np.fromstring(wave_piece.findtext(
            "./Polys/DataArray[@Name='offsets']"
        ), sep=" ", dtype=int)
        if (points.shape != (2000, 3) or
                wave_piece.get("NumberOfPoints") != "2000" or
                wave_piece.get("NumberOfPolys") != "999" or
                wave_cells.shape != (999 * 4,) or
                wave_offsets.shape != (999,) or
                np.any(wave_cells < 0) or np.any(wave_cells >= 2000) or
                not np.array_equal(wave_offsets, 4 * np.arange(1, 1000))):
            raise ValueError(f"wave frame {number} has wrong geometry")
        if frame == 0:
            wave_connectivity = wave_cells
            saved_wave_offsets = wave_offsets
        elif not np.array_equal(wave_connectivity, wave_cells):
            raise ValueError(f"wave frame {number} topology changed")
        wave_samples[frame] = points[WAVE_SAMPLE_INDICES]
        wave_extrema[frame] = (np.min(points[:, 2]), np.max(points[:, 2]))
        for body, name in enumerate(("float", "spar")):
            source = ElementTree.parse(
                scene / f"body{body + 1}_{name}" / f"{name}_{number}.vtp"
            ).getroot()
            piece = source.find("./PolyData/Piece")
            points = _points(source, "./PolyData/Piece/Points/DataArray")
            connectivity = np.fromstring(piece.findtext(
                "./Polys/DataArray[@Name='connectivity']"
            ), sep=" ", dtype=int)
            offsets = np.fromstring(piece.findtext(
                "./Polys/DataArray[@Name='offsets']"
            ), sep=" ", dtype=int)
            cell_area = np.fromstring(piece.findtext(
                "./CellData/DataArray[@Name='CellArea']"
            ), sep=" ")
            facets = (516, 356)[body]
            if (points.shape != ((258, 180)[body], 3) or
                    piece.get("NumberOfPoints") != str((258, 180)[body]) or
                    piece.get("NumberOfPolys") != str(facets) or
                    connectivity.shape != (facets * 3,) or
                    offsets.shape != (facets,) or
                    cell_area.shape != (facets,) or
                    np.any(connectivity < 0) or
                    np.any(connectivity >= len(points)) or
                    not np.array_equal(offsets, 3 * np.arange(1, facets + 1)) or
                    not np.isfinite(cell_area).all() or
                    np.any(cell_area <= 0)):
                raise ValueError(f"body {body + 1} frame {number} changed")
            triangles = points[connectivity.reshape(-1, 3)]
            centers = triangles.mean(axis=1)
            area = np.linalg.norm(np.cross(
                triangles[:, 1] - triangles[:, 0],
                triangles[:, 2] - triangles[:, 0],
            ), axis=1) / 2
            body_stats[body, frame] = np.r_[
                centers.mean(axis=0), points.min(axis=0),
                points.max(axis=0), area.sum(), cell_area.sum(),
            ]
        mooring = ElementTree.parse(
            scene / "mooring1" / f"mooring_{number}.vtp"
        ).getroot()
        pieces = mooring.findall("./PolyData/Piece")
        if len(pieces) != 3:
            raise ValueError(f"mooring frame {number} has wrong line count")
        for line, piece in enumerate(pieces):
            points = _points(piece, "./Points/DataArray")
            connectivity = np.fromstring(piece.findtext(
                "./Lines/DataArray[@Name='connectivity']"
            ), sep=" ", dtype=int).reshape(-1, 2)
            offsets = np.fromstring(piece.findtext(
                "./Lines/DataArray[@Name='offsets']"
            ), sep=" ", dtype=int)
            tension = np.fromstring(piece.findtext(
                "./CellData/DataArray[@Name='Segment Tension']"
            ), sep=" ")
            expected_connectivity = np.column_stack((
                np.arange(20), np.arange(1, 21),
            ))
            if (piece.get("NumberOfPoints") != "21" or
                    piece.get("NumberOfLines") != "20" or
                    points.shape != (21, 3) or tension.shape != (20,)
                    or not np.isfinite(tension).all()
                    or not np.array_equal(connectivity, expected_connectivity)
                    or not np.array_equal(offsets, 2 * np.arange(1, 21))):
                raise ValueError(f"mooring frame {number} line {line} changed")
            line_nodes[line, frame] = points
            line_tension[line, frame] = tension
    return {
        "scene_wave_sample_indices": WAVE_SAMPLE_INDICES,
        "scene_wave_connectivity": wave_connectivity,
        "scene_wave_offsets": saved_wave_offsets,
        "scene_wave_samples": wave_samples,
        "scene_wave_z_extrema": wave_extrema,
        "scene_body_stats": body_stats,
        "scene_line_nodes": line_nodes,
        "scene_line_tension": line_tension,
    }


def main() -> None:
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--full-fine-dir", type=Path, required=True)
    parser.add_argument("--corrected-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    arms = [
        _source(args.full_fine_dir, "fullfine_000625_", stride=10),
        _source(args.full_fine_dir, "fullfine_0003125_", stride=10),
        _source(args.corrected_dir),
    ]
    for arm in arms:
        grid = np.arange(801) * .1
        for kind in ("wave", "body1", "body2", "pto", "mooring",
                     "fairlead", "line1", "line2", "line3"):
            np.testing.assert_allclose(arm[kind][:, 0], grid, rtol=0,
                                       atol=1e-8)
        np.testing.assert_allclose(arm["fairlead_dense"][:, 0],
                                   np.arange(8001) * .01, rtol=0, atol=1e-8)
    result = {
        "time": np.arange(801) * .1,
        "phase": np.stack([arm["phase"][:, 0] for arm in arms]),
        "wave": np.stack([arm["wave"][:, 1] for arm in arms]),
        "body_position": np.stack([
            np.stack([arm["body1"][:, 1:7], arm["body2"][:, 1:7]])
            for arm in arms
        ]),
        "body_velocity": np.stack([
            np.stack([arm["body1"][:, 7:13], arm["body2"][:, 7:13]])
            for arm in arms
        ]),
        "pto_force": np.stack([arm["pto"][:, 15] for arm in arms]),
        "pto_power": np.stack([arm["pto"][:, 21] for arm in arms]),
        "connection_position": np.stack([
            arm["mooring"][:, 1:7] for arm in arms
        ]),
        "connection_velocity": np.stack([
            arm["mooring"][:, 7:13] for arm in arms
        ]),
        "connection_force": np.stack([
            arm["mooring"][:, 13:19] for arm in arms
        ]),
        "fairlead_tension": np.stack([
            arm["fairlead"][:, 1:4] for arm in arms
        ]),
        "fairlead_dense_time": np.arange(8001) * .01,
        "fairlead_dense_tension": np.stack([
            arm["fairlead_dense"][:, 1:4] for arm in arms
        ]),
        "line_nodes": np.stack([
            np.stack([
                arm[f"line{line}"][:, 1:64].reshape(801, 21, 3)
                for line in (1, 2, 3)
            ]) for arm in arms
        ]),
        "line_tension": np.stack([
            np.stack([
                arm[f"line{line}"][:, 64:84]
                for line in (1, 2, 3)
            ]) for arm in arms
        ]),
    }
    result.update(_scene(args.corrected_dir / "scene"))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.output, **result)
    print(f"saved {args.output} ({args.output.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()

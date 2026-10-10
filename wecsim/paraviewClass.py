"""Wave-surface and body-mesh VTP output for ParaView.

The upstream file contained unfinished MATLAB syntax and could not be imported.
"""

from pathlib import Path
from xml.etree import ElementTree

import numpy as np


class ParaviewClass:
    def __init__(self, waves):
        self.waves = waves

    def waveElevationGrid(self, t, X, Y):
        """Return the undisturbed wave elevation over an X/Y grid at time t.

        This follows ``waveClass.waveElevationGrid`` in MATLAB WEC-Sim for
        no-wave, regular, and frequency-resolved irregular waves. The arrays
        X and Y must have matching or broadcastable shapes.
        """
        X, Y = np.broadcast_arrays(np.asarray(X, dtype=float),
                                   np.asarray(Y, dtype=float))
        waves = self.waves
        if waves.wType in ("noWave", "noWaveCIC"):
            return np.zeros_like(X)

        if waves.wType in ("regular", "regularCIC"):
            direction = np.deg2rad(np.asarray(waves.waveDir).item())
            position = X * np.cos(direction) + Y * np.sin(direction)
            return np.asarray(waves.A).item() * np.cos(
                -np.asarray(waves.k).item() * position
                + np.asarray(waves.w).item() * t
            )

        if waves.wType in ("irregular", "spectrumImport"):
            frequencies = np.atleast_1d(waves.w)
            wavenumbers = np.atleast_1d(waves.k)
            amplitudes = np.atleast_1d(waves.A)
            intervals = np.atleast_1d(waves.dw)
            directions = np.atleast_1d(waves.waveDir)
            spreads = np.atleast_1d(waves.waveSpread)
            phases = np.asarray(waves.phase)
            if phases.shape != (len(directions), len(frequencies)):
                raise ValueError("phase must have one row per direction and one column per frequency")
            if not (len(frequencies) == len(wavenumbers) == len(amplitudes) == len(intervals)):
                raise ValueError("frequency, wavenumber, amplitude, and interval lengths must match")
            if len(directions) != len(spreads):
                raise ValueError("direction and spread lengths must match")

            elevation = np.zeros_like(X)
            for direction_index, direction in enumerate(directions):
                angle = np.deg2rad(direction)
                position = X * np.cos(angle) + Y * np.sin(angle)
                for frequency_index, frequency in enumerate(frequencies):
                    amplitude = np.sqrt(
                        amplitudes[frequency_index]
                        * intervals[frequency_index]
                        * spreads[direction_index]
                    )
                    elevation += amplitude * np.cos(
                        -wavenumbers[frequency_index] * position
                        + frequency * t
                        + phases[direction_index, frequency_index]
                    )
            return elevation

        raise NotImplementedError(
            f"wave-surface visualization is not implemented for {waves.wType!r}"
        )

    def write_paraview_vtp_wave(
        self, times, directory, *, domain_size,
        num_points_x=None, num_points_y=None, water_depth=None,
        num_moordyn=0,
    ) -> tuple[Path, ...]:
        """Write a WEC-Sim-style wave grid as numbered VTP PolyData frames.

        The grid covers ``[-domain_size, domain_size]`` on both axes. This
        writes the wave and ``ground.txt`` only; body and mooring meshes are
        separate visualization outputs.
        """
        times = np.asarray(times, dtype=float)
        if (times.ndim != 1 or not len(times)
                or not np.isfinite(times).all()
                or np.any(np.diff(times) <= 0)):
            raise ValueError("wave VTP times must be finite and increasing")
        nx = self.waves.viz["numPointsX"] if num_points_x is None else num_points_x
        ny = self.waves.viz["numPointsY"] if num_points_y is None else num_points_y
        if (isinstance(nx, bool) or not isinstance(nx, (int, np.integer))
                or nx < 2 or isinstance(ny, bool)
                or not isinstance(ny, (int, np.integer)) or ny < 2):
            raise ValueError("wave VTP grid needs at least two points per axis")
        if (not np.isfinite(domain_size) or domain_size <= 0
                or isinstance(num_moordyn, bool)
                or not isinstance(num_moordyn, (int, np.integer))
                or num_moordyn < 0):
            raise ValueError("wave VTP domain and mooring count must be valid")
        depth = self.waves.waterDepth if water_depth is None else water_depth
        depth = np.asarray(depth, dtype=float)
        if depth.size != 1 or not np.isfinite(depth).all() or depth.item() <= 0:
            raise ValueError("wave VTP needs a positive water depth")

        directory = Path(directory)
        wave_dir = directory / "waves"
        wave_dir.mkdir(parents=True, exist_ok=True)
        (directory / "ground.txt").write_text(
            f"{domain_size:g}\n{depth.item():g}\n{num_moordyn}\n",
            encoding="utf-8",
        )
        x = np.linspace(-domain_size, domain_size, nx)
        y = np.linspace(-domain_size, domain_size, ny)
        X, Y = np.meshgrid(x, y)
        connectivity = np.array([
            (row * nx + col, row * nx + col + 1,
             (row + 1) * nx + col + 1, (row + 1) * nx + col)
            for row in range(ny - 1) for col in range(nx - 1)
        ])
        paths = []
        for index, time in enumerate(times, start=1):
            Z = self.waveElevationGrid(float(time), X, Y)
            if Z.shape != X.shape or not np.isfinite(Z).all():
                raise ValueError("wave VTP elevation must be finite on the grid")
            root = ElementTree.Element("VTKFile", type="PolyData", version="0.1")
            polydata = ElementTree.SubElement(root, "PolyData")
            piece = ElementTree.SubElement(
                polydata, "Piece", NumberOfPoints=str(nx * ny),
                NumberOfPolys=str((nx - 1) * (ny - 1)),
            )
            points = ElementTree.SubElement(piece, "Points")
            point_values = ElementTree.SubElement(
                points, "DataArray", type="Float32",
                NumberOfComponents="3", format="ascii",
            )
            point_values.text = "\n" + "\n".join(
                f"{xx:.5f} {yy:.5f} {zz:.5f}"
                for xx, yy, zz in zip(X.ravel(), Y.ravel(), Z.ravel())
            ) + "\n"
            polys = ElementTree.SubElement(piece, "Polys")
            cells = ElementTree.SubElement(
                polys, "DataArray", type="Int32", Name="connectivity",
                format="ascii",
            )
            cells.text = "\n" + "\n".join(
                " ".join(map(str, cell)) for cell in connectivity
            ) + "\n"
            offsets = ElementTree.SubElement(
                polys, "DataArray", type="Int32", Name="offsets",
                format="ascii",
            )
            offsets.text = " " + " ".join(
                str(4 * cell) for cell in range(1, len(connectivity) + 1)
            ) + " "
            path = wave_dir / f"waves_{index}.vtp"
            ElementTree.ElementTree(root).write(
                path, encoding="utf-8", xml_declaration=True,
            )
            paths.append(path)
        return tuple(paths)

    def write_paraview_vtp_mooring(
        self, source_times, frame_times, directory, *, lines, tensions,
        mooring_index=1,
    ) -> tuple[Path, ...]:
        """Write MoorDyn node paths and segment tensions as numbered VTP frames.

        Each line has shape ``(source_time, node, xyz)`` and its tensions have
        shape ``(source_time, node - 1)``. Values are linearly interpolated to
        ``frame_times`` as in the pinned MATLAB ``writeParaviewMooring``.
        """
        source_times = np.asarray(source_times, dtype=float)
        frame_times = np.asarray(frame_times, dtype=float)
        if (source_times.ndim != 1 or len(source_times) < 2
                or not np.isfinite(source_times).all()
                or np.any(np.diff(source_times) <= 0)
                or frame_times.ndim != 1 or not len(frame_times)
                or not np.isfinite(frame_times).all()
                or np.any(np.diff(frame_times) <= 0)
                or frame_times[0] < source_times[0]
                or frame_times[-1] > source_times[-1]
                or isinstance(mooring_index, bool)
                or not isinstance(mooring_index, (int, np.integer))
                or mooring_index < 1):
            raise ValueError("mooring VTP needs ordered in-range times and a positive index")
        if len(lines) == 0 or len(lines) != len(tensions):
            raise ValueError("mooring VTP needs node and tension data for each line")

        interpolated = []
        for points, load in zip(lines, tensions):
            points = np.asarray(points, dtype=float)
            load = np.asarray(load, dtype=float)
            if (points.ndim != 3 or points.shape[0] != len(source_times)
                    or points.shape[1] < 2 or points.shape[2] != 3
                    or load.shape != (len(source_times), points.shape[1] - 1)
                    or not np.isfinite(points).all()
                    or not np.isfinite(load).all()):
                raise ValueError("mooring VTP lines need finite nodes and segment tensions")
            moved = np.empty((len(frame_times), points.shape[1], 3))
            forces = np.empty((len(frame_times), points.shape[1] - 1))
            for node in range(points.shape[1]):
                for axis in range(3):
                    moved[:, node, axis] = np.interp(
                        frame_times, source_times, points[:, node, axis],
                    )
            for segment in range(points.shape[1] - 1):
                forces[:, segment] = np.interp(
                    frame_times, source_times, load[:, segment],
                )
            interpolated.append((moved, forces))

        line_dir = Path(directory) / f"mooring{mooring_index}"
        line_dir.mkdir(parents=True, exist_ok=True)
        paths = []
        for frame in range(len(frame_times)):
            root = ElementTree.Element("VTKFile", type="PolyData", version="0.1")
            polydata = ElementTree.SubElement(root, "PolyData")
            for nodes, tension in interpolated:
                points = nodes[frame]
                forces = tension[frame]
                piece = ElementTree.SubElement(
                    polydata, "Piece", NumberOfPoints=str(len(points)),
                    NumberOfLines=str(len(forces)),
                )
                point_data = ElementTree.SubElement(
                    ElementTree.SubElement(piece, "Points"), "DataArray",
                    type="Float32", NumberOfComponents="3", format="ascii",
                )
                point_data.text = "\n" + "\n".join(
                    f"{x:.5f} {y:.5f} {z:.5f}" for x, y, z in points
                ) + "\n"
                cells = ElementTree.SubElement(piece, "Lines")
                connectivity = ElementTree.SubElement(
                    cells, "DataArray", type="Int32", Name="connectivity",
                    format="ascii",
                )
                connectivity.text = "\n" + "\n".join(
                    f"{segment} {segment + 1}"
                    for segment in range(len(forces))
                ) + "\n"
                offsets = ElementTree.SubElement(
                    cells, "DataArray", type="Int32", Name="offsets",
                    format="ascii",
                )
                offsets.text = " " + " ".join(
                    str(2 * (segment + 1)) for segment in range(len(forces))
                ) + " "
                cell_data = ElementTree.SubElement(piece, "CellData")
                tension_data = ElementTree.SubElement(
                    cell_data, "DataArray", type="Float32",
                    Name="Segment Tension", NumberOfComponents="1",
                    format="ascii",
                )
                tension_data.text = " " + " ".join(
                    f"{value:.8g}" for value in forces
                ) + " "
            path = line_dir / f"mooring_{frame + 1}.vtp"
            ElementTree.ElementTree(root).write(
                path, encoding="utf-8", xml_declaration=True,
            )
            paths.append(path)
        return tuple(paths)

    def write_paraview_vtp(
        self, times, directory, *, body_name, vertices, faces, poses,
        body_index=1, cell_areas=None, hydrostatic_pressure=None,
        nonlinear_wave_pressure=None, linear_wave_pressure=None,
    ) -> tuple[Path, ...]:
        """Write numbered body VTP meshes using WEC-Sim's XYZ pose order.

        ``faces`` use zero-based triangular vertex indices. Each row of
        ``poses`` holds surge, sway, heave, roll, pitch, and yaw in SI units.
        Optional pressure arrays have one row per time and one column per face.
        """
        times = np.asarray(times, dtype=float)
        vertices = np.asarray(vertices, dtype=float)
        faces = np.asarray(faces)
        poses = np.asarray(poses, dtype=float)
        if (times.ndim != 1 or not len(times)
                or not np.isfinite(times).all()
                or np.any(np.diff(times) <= 0)
                or vertices.ndim != 2 or vertices.shape[1] != 3
                or not len(vertices) or not np.isfinite(vertices).all()
                or faces.ndim != 2 or faces.shape[1] != 3
                or not len(faces) or not np.issubdtype(faces.dtype, np.integer)
                or np.any(faces < 0) or np.any(faces >= len(vertices))
                or poses.shape != (len(times), 6)
                or not np.isfinite(poses).all()):
            raise ValueError("body VTP needs finite time, vertices, triangles, and six-DOF poses")
        if (not isinstance(body_name, str) or not body_name
                or body_name in (".", "..") or "/" in body_name
                or "\\" in body_name or isinstance(body_index, bool)
                or not isinstance(body_index, (int, np.integer))
                or body_index < 1):
            raise ValueError("body VTP needs a name and positive body index")
        if cell_areas is None:
            edges = vertices[faces[:, 1]] - vertices[faces[:, 0]]
            sides = vertices[faces[:, 2]] - vertices[faces[:, 0]]
            areas = np.linalg.norm(np.cross(edges, sides), axis=1) / 2
        else:
            areas = np.asarray(cell_areas, dtype=float)
        if (areas.shape != (len(faces),) or not np.isfinite(areas).all()
                or np.any(areas < 0)):
            raise ValueError("body VTP cell areas must be finite and nonnegative")
        pressure_fields = {}
        for name, values in (
            ("HydrostaticPressure", hydrostatic_pressure),
            ("WavePressureNonLinear", nonlinear_wave_pressure),
            ("WavePressureLinear", linear_wave_pressure),
        ):
            if values is not None:
                field = np.asarray(values, dtype=float)
                if (field.shape != (len(times), len(faces))
                        or not np.isfinite(field).all()):
                    raise ValueError(f"body VTP {name} must match time and faces")
                pressure_fields[name] = field

        body_dir = Path(directory) / f"body{body_index}_{body_name}"
        body_dir.mkdir(parents=True, exist_ok=True)
        paths = []
        for index, pose in enumerate(poses, start=1):
            roll, pitch, yaw = pose[3:]
            cx, sx = np.cos(roll), np.sin(roll)
            cy, sy = np.cos(pitch), np.sin(pitch)
            cz, sz = np.cos(yaw), np.sin(yaw)
            rx = np.array([[1, 0, 0], [0, cx, -sx], [0, sx, cx]])
            ry = np.array([[cy, 0, sy], [0, 1, 0], [-sy, 0, cy]])
            rz = np.array([[cz, -sz, 0], [sz, cz, 0], [0, 0, 1]])
            moved = vertices @ (rz @ ry @ rx).T + pose[:3]
            root = ElementTree.Element("VTKFile", type="PolyData", version="0.1")
            polydata = ElementTree.SubElement(root, "PolyData")
            piece = ElementTree.SubElement(
                polydata, "Piece", NumberOfPoints=str(len(vertices)),
                NumberOfPolys=str(len(faces)),
            )
            points = ElementTree.SubElement(piece, "Points")
            point_values = ElementTree.SubElement(
                points, "DataArray", type="Float32",
                NumberOfComponents="3", format="ascii",
            )
            point_values.text = "\n" + "\n".join(
                f"{x:.5f} {y:.5f} {z:.5f}" for x, y, z in moved
            ) + "\n"
            polys = ElementTree.SubElement(piece, "Polys")
            cells = ElementTree.SubElement(
                polys, "DataArray", type="Int32", Name="connectivity",
                format="ascii",
            )
            cells.text = "\n" + "\n".join(
                " ".join(map(str, face)) for face in faces
            ) + "\n"
            offsets = ElementTree.SubElement(
                polys, "DataArray", type="Int32", Name="offsets",
                format="ascii",
            )
            offsets.text = " " + " ".join(
                str(3 * cell) for cell in range(1, len(faces) + 1)
            ) + " "
            cell_data = ElementTree.SubElement(piece, "CellData")
            for name, values in (("CellArea", areas), *pressure_fields.items()):
                field = values if name == "CellArea" else values[index - 1]
                data = ElementTree.SubElement(
                    cell_data, "DataArray", type="Float32", Name=name,
                    NumberOfComponents="1", format="ascii",
                )
                data.text = " " + " ".join(f"{value:.8g}" for value in field) + " "
            path = body_dir / f"{body_name}_{index}.vtp"
            ElementTree.ElementTree(root).write(
                path, encoding="utf-8", xml_declaration=True,
            )
            paths.append(path)
        return tuple(paths)

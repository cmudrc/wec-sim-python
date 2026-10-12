"""Build constant small-motion body maps from named device coordinates."""

import re
from typing import Mapping

import numpy as np


AXES = ("surge", "sway", "heave", "roll", "pitch", "yaw")


def _pivot(value, center):
    if not isinstance(value, Mapping) or len(value) != 1:
        raise ValueError("rotation pivot needs world or body-local point")
    if "world" in value:
        raw = value["world"]
        origin = np.zeros(3)
    elif "point" in value:
        raw = value["point"]
        origin = center
    else:
        raise ValueError("rotation pivot needs world or body-local point")
    point = np.asarray(raw, dtype=float)
    if point.shape != (3,) or not np.isfinite(point).all():
        raise ValueError("rotation pivot must be a finite three-vector")
    return origin + point


def build_coordinate_maps(constraint, bodies, body_names, centers):
    """Return six-by-N body maps and names of their independent coordinates.

    A coordinate can move one or several named bodies. Reusing one coordinate
    for multiple body DOFs expresses a shared linearized motion. The solver
    checks the assembled mass matrix for singular or redundant coordinates.
    """
    if len(centers) != len(bodies):
        raise ValueError("body centers must align with body definitions")
    if "coordinates" not in constraint:
        maps = []
        for body_spec in bodies:
            mapping = np.asarray(body_spec.get("coordinate_map"), dtype=float)
            if (mapping.ndim != 2 or mapping.shape[0] != 6
                    or mapping.shape[1] < 1 or not np.isfinite(mapping).all()):
                raise ValueError("body.coordinate_map must be a finite 6-by-N matrix")
            maps.append(mapping)
        count = maps[0].shape[1]
        if any(mapping.shape[1] != count for mapping in maps):
            raise ValueError("all coordinate maps must use the same coordinate count")
        return maps, tuple(f"coordinate{index}" for index in range(1, count + 1))

    if any("coordinate_map" in body for body in bodies):
        raise ValueError("use named coordinates or body.coordinate_map, not both")
    definitions = constraint["coordinates"]
    if not isinstance(definitions, list) or not definitions:
        raise ValueError("constraint.coordinates must be a nonempty list")
    maps = [np.zeros((6, len(definitions))) for _ in bodies]
    names = []
    for coordinate, definition in enumerate(definitions):
        if not isinstance(definition, Mapping):
            raise ValueError("each coordinate must be an object")
        if set(definition) != {"name", "motions"}:
            raise ValueError("each coordinate needs only name and motions")
        name = definition["name"]
        if (not isinstance(name, str)
                or re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", name) is None
                or name in names):
            raise ValueError("coordinate names must be unique identifiers")
        names.append(name)
        motions = definition["motions"]
        if not isinstance(motions, list) or not motions:
            raise ValueError("each coordinate needs at least one body motion")
        seen = set()
        for motion in motions:
            if (not isinstance(motion, Mapping) or {"body", "dof"} - motion.keys()
                    or motion.keys() - {"body", "dof", "scale", "pivot"}):
                raise ValueError("coordinate motion needs body, dof, optional scale and pivot")
            body = motion["body"]
            if isinstance(body, bool):
                raise ValueError("coordinate motion body is unknown")
            if isinstance(body, int):
                body_index = body - 1
            elif isinstance(body, str) and body in body_names:
                body_index = body_names.index(body)
            else:
                raise ValueError("coordinate motion body is unknown")
            if body_index < 0 or body_index >= len(bodies):
                raise ValueError("coordinate motion body is unknown")
            if motion["dof"] not in AXES:
                raise ValueError("coordinate motion dof is unsupported")
            dof = AXES.index(motion["dof"])
            if (body_index, dof) in seen:
                raise ValueError("a coordinate cannot repeat one body motion")
            seen.add((body_index, dof))
            scale = motion.get("scale", 1)
            if isinstance(scale, bool):
                raise ValueError("coordinate motion scale must be finite and nonzero")
            try:
                scale = float(scale)
            except (TypeError, ValueError) as exc:
                raise ValueError("coordinate motion scale must be finite and nonzero") from exc
            if not np.isfinite(scale) or scale == 0:
                raise ValueError("coordinate motion scale must be finite and nonzero")
            maps[body_index][dof, coordinate] += scale
            if "pivot" in motion:
                if dof < 3:
                    raise ValueError("only rotational motions use a pivot")
                pivot = _pivot(motion["pivot"], np.asarray(centers[body_index]))
                axis = np.eye(3)[dof - 3]
                maps[body_index][:3, coordinate] += scale * np.cross(
                    axis, np.asarray(centers[body_index]) - pivot,
                )
    return maps, tuple(names)


def initial_coordinate(value, names, label):
    """Accept a vector or a sparse dictionary keyed by coordinate name."""
    if isinstance(value, Mapping):
        unknown = value.keys() - set(names)
        if unknown:
            raise ValueError(f"{label} has unknown coordinates {sorted(unknown)}")
        values = [value.get(name, 0) for name in names]
    else:
        values = value
    coordinates = np.asarray(values, dtype=float)
    if (coordinates.shape != (len(names),)
            or not np.isfinite(coordinates).all()):
        raise ValueError(f"{label} must have one finite value per coordinate")
    return coordinates

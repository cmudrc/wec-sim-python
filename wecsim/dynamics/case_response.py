"""Result arrays returned by a case simulation."""

from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class CaseResponse:
    time: np.ndarray
    body_position: np.ndarray  # (time, body, six WEC-Sim coordinates)
    body_velocity: np.ndarray
    hydro_files: tuple[Path | None, ...]
    pto_force: np.ndarray | None = None
    pto_label: str | None = None
    wave_elevation: np.ndarray | None = None
    total_heave_force: np.ndarray | None = None
    auxiliary_files: tuple[Path, ...] = ()
    pto_generalized_force: np.ndarray | None = None
    extra_outputs: tuple[tuple[str, np.ndarray], ...] = ()
    coordinate_position: np.ndarray | None = None
    coordinate_velocity: np.ndarray | None = None
    pto_stroke: np.ndarray | None = None
    pto_velocity: np.ndarray | None = None
    pto_absorbed_power: np.ndarray | None = None

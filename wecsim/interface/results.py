"""Public motion, PTO, and WEC result records."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from wecsim.dynamics.caseDynamics import CaseResponse


@dataclass(frozen=True)
class MotionHistory:
    position: np.ndarray
    velocity: np.ndarray


@dataclass(frozen=True)
class FlexibleModeHistory:
    position: np.ndarray
    velocity: np.ndarray
    acceleration: np.ndarray


@dataclass(frozen=True)
class DirectDriveHistory:
    shaft_velocity: np.ndarray
    shaft_torque: np.ndarray
    inertia_torque: np.ndarray
    friction_torque: np.ndarray
    generator_torque: np.ndarray
    current: np.ndarray
    voltage: np.ndarray
    resistance_loss: np.ndarray
    electrical_power: np.ndarray
    mechanical_power: np.ndarray


@dataclass(frozen=True)
class LinearGeneratorHistory:
    flux_d: np.ndarray
    flux_q: np.ndarray
    electrical_angle: np.ndarray
    friction_force: np.ndarray
    electrical_power: np.ndarray
    phase_current: np.ndarray
    phase_voltage: np.ndarray


@dataclass(frozen=True)
class PTOHistory:
    """PTO stroke in metres, or angle in radians for a rotational PTO.

    Force and absorbed power include a configured velocity-control gain;
    the optional control histories expose that contribution separately.
    """

    stroke: np.ndarray
    velocity: np.ndarray
    force: np.ndarray
    absorbed_power: np.ndarray
    direct_drive: DirectDriveHistory | None = None
    linear_generator: LinearGeneratorHistory | None = None
    control_force: np.ndarray | None = None
    control_absorbed_power: np.ndarray | None = None


@dataclass(frozen=True)
class WECResult:
    time: np.ndarray
    bodies: dict[str, MotionHistory]
    coordinates: dict[str, MotionHistory]
    ptos: dict[str, PTOHistory]
    wave_elevation: np.ndarray | None
    case: dict
    raw: CaseResponse
    flexible_modes: dict[str, FlexibleModeHistory] = field(default_factory=dict)
    body_forces: dict[str, np.ndarray] = field(default_factory=dict)

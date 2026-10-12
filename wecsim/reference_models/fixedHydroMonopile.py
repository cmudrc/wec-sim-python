"""Stationary hydro and fixed-joint outputs for the published monopile case."""

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from wecsim.hydrodynamics.bodyClass import BodyClass
from wecsim.environment.irregularWave import IrregularComponents, synthesize_irregular_response


@dataclass(frozen=True)
class FixedHydroMonopileResponse:
    time: np.ndarray
    wave_elevation: np.ndarray
    excitation_force: np.ndarray
    body_position: np.ndarray
    body_velocity: np.ndarray
    body_force_total: np.ndarray
    joint_force: np.ndarray


def run_fixed_hydro_monopile(
    hydro_file: str | Path,
    components: IrregularComponents,
    *,
    tower_mass: float = 1_031_930.,
    tower_center: tuple[float, float, float] = (0., 0., 25.),
    joint_center: tuple[float, float, float] = (0., 0., 0.),
    dt: float = .01,
    end_time: float = 400.,
    ramp_time: float = 100.,
    rho: float = 1025.,
    g: float = 9.81,
) -> FixedHydroMonopileResponse:
    """Evaluate the published fixed hydro body and attached tower.

    The monopile hydro body's equilibrium weight and buoyancy cancel. Both
    bodies are stationary, so radiation, added-mass, drag, and acceleration
    forces vanish. ``joint_force`` follows the published model's two logged
    fixed-joint reaction conventions; it is not a general support solver.
    """
    if not isinstance(components, IrregularComponents):
        raise TypeError("monopile needs irregular wave components")
    scalars = np.asarray((tower_mass, dt, end_time, ramp_time, rho, g),
                         dtype=float)
    tower = np.asarray(tower_center, dtype=float)
    joint = np.asarray(joint_center, dtype=float)
    if (not np.isfinite(scalars).all() or tower_mass <= 0 or dt <= 0
            or end_time <= 0 or ramp_time < 0 or rho <= 0 or g <= 0
            or tower.shape != (3,) or joint.shape != (3,)
            or not np.isfinite(tower).all() or not np.isfinite(joint).all()):
        raise ValueError("monopile body, joint, and simulation settings are invalid")
    steps = round(end_time / dt)
    if not np.isclose(steps * dt, end_time, rtol=0, atol=1e-10):
        raise ValueError("monopile duration must have an integer number of steps")

    body = BodyClass(str(hydro_file))
    body.bodyNumber = 1
    body.readH5file()
    if int(np.asarray(body.dof).item()) != 6:
        raise ValueError("monopile hydrodynamics need six rigid DOFs")
    hydro_center = np.asarray(body.cg, dtype=float).ravel()
    sea = synthesize_irregular_response(
        hydro_file, components, dt=dt, end_time=end_time,
        ramp_time=ramp_time, body_number=1, rho=rho, g=g,
    )
    force = sea.excitation_force[:, :6]
    count = len(sea.time)
    position = np.zeros((count, 2, 6))
    position[:, 0, :3] = hydro_center
    position[:, 1, :3] = tower
    velocity = np.zeros_like(position)
    body_force = np.zeros_like(position)
    body_force[:, 0] = force
    body_force[:, 1, 2] = -tower_mass * g

    joint_force = np.zeros_like(position)
    joint_force[:, 0, :3] = -force[:, :3] - body_force[:, 1, :3]
    joint_force[:, 0, 3:] = (
        -force[:, 3:]
        - np.cross(hydro_center - joint, force[:, :3])
        - np.cross(tower - joint, body_force[:, 1, :3])
    )
    joint_force[:, 1] = -body_force[:, 1]
    return FixedHydroMonopileResponse(
        sea.time, sea.elevation, force, position, velocity,
        body_force, joint_force,
    )

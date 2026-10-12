"""Three-coordinate WaveBot motion for the published impedance experiment."""

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.io import loadmat

from wecsim.hydrodynamics.bodyClass import BodyClass
from wecsim.dynamics.generalDynamics import BodyMotion, DynamicBody, GeneralizedDynamics


@dataclass(frozen=True)
class WaveBotImpedanceResponse:
    time: np.ndarray
    body_position: np.ndarray
    body_velocity: np.ndarray
    body_acceleration: np.ndarray
    applied_force: np.ndarray
    mooring_force: np.ndarray
    linear_damping_force: np.ndarray
    drag_force: np.ndarray


def run_wavebot_impedance(
    hydro_file: str | Path,
    multisine_file: str | Path,
    *,
    dt: float = .01,
    end_time: float = 600.,
    rho: float = 1025.,
    g: float = 9.81,
    source_linear_damping: bool = False,
) -> WaveBotImpedanceResponse:
    """Run the published no-wave, multisine-forced WaveBot experiment.

    The rigid body moves in surge, heave, and pitch. Radiation uses the
    WAMIT-derived 10 s impulse response; added mass remains implicit. The
    mooring, quadratic drag, and multisine gains match the pinned
    ``CalcImpedance`` input. By default, the intended surge/heave/pitch linear
    damping acts on each matching velocity. The published MATLAB input uses
    linear indexing and puts all three coefficients in the surge-velocity
    column; select ``source_linear_damping`` to reproduce that input exactly.
    This runner does not implement the later adaptive controller.
    """
    if (not np.isfinite([dt, end_time, rho, g]).all() or dt <= 0
            or end_time <= 0 or rho <= 0 or g <= 0):
        raise ValueError("WaveBot simulation settings must be finite and positive")
    if not isinstance(source_linear_damping, bool):
        raise TypeError("source_linear_damping must be a boolean")
    steps = round(end_time / dt)
    if not np.isclose(steps * dt, end_time, rtol=0, atol=1e-10):
        raise ValueError("WaveBot duration must contain an integer number of steps")
    time = np.arange(steps + 1) * dt
    mat = loadmat(multisine_file, variable_names=["multisine"])
    if "multisine" not in mat:
        raise ValueError("multisine_file must contain the published multisine struct")
    record = mat["multisine"][0, 0]
    signal_time = np.asarray(record["time"], dtype=float).ravel()
    signal_values = np.asarray(record["signals"][0, 0]["values"], dtype=float)
    if (signal_values.ndim != 2 or signal_values.shape != (len(signal_time), 3)
            or len(signal_time) < 2 or signal_time[0] != 0
            or signal_time[-1] < end_time - .01
            or not np.all(np.diff(signal_time) > 0)
            or not np.isfinite(signal_time).all()
            or not np.isfinite(signal_values).all()):
        raise ValueError("multisine must have increasing time and three finite channels")
    applied = np.column_stack([
        np.interp(time, signal_time, signal_values[:, axis])
        for axis in range(3)
    ])
    late = time > signal_time[-1]
    applied[late] = (
        signal_values[-1]
        + (time[late, None] - signal_time[-1])
        * (signal_values[-1] - signal_values[-2])
        / (signal_time[-1] - signal_time[-2])
    )
    applied *= np.array([50., 100., 4.])
    # The published userDefinedFunctions.m negates the third logged command
    # before using it as the pitch torque in impedance identification.
    applied[:, 2] *= -1

    body = BodyClass(str(hydro_file))
    body.bodyNumber = 1
    body.bodyTotal = 1
    body.readH5file()
    if int(np.asarray(body.dof).item()) != 6:
        raise ValueError("WaveBot HDF5 must contain one six-DOF body")
    body.mass = 1156.5
    zero = np.zeros((6, 6))
    body.hydroStiffness = zero.copy()
    body.viscDrag = {
        "Drag": zero.copy(), "cd": np.zeros(6),
        "characteristicArea": np.zeros(6),
    }
    body.linearDamping = zero.copy()
    memory_time = np.arange(round(10 / dt) + 1) * dt
    body.hydroForcePre(
        [], [0], len(memory_time), memory_time, [], dt, rho, g,
        "noWaveCIC", np.vstack((time, np.zeros_like(time))),
        1, 1, 0, 0, 0,
    )
    hydro = body.hydroForce
    center = np.asarray(body.cg, dtype=float).ravel()
    mapping = np.zeros((6, 3))
    mapping[0, 0] = mapping[2, 1] = mapping[4, 2] = 1
    inertia = np.diag([1156.5] * 3 + [84.] * 3)
    linear = np.diag([1000., 0., 1000., 0., 100., 0.])
    if source_linear_damping:
        linear = np.zeros((6, 6))
        linear[[0, 2, 4], 0] = [1000., 1000., 100.]
    drag = .5 * rho * np.array([1.15, 1.15, 1., .5, .5, 0.]) * np.array(
        [2.9568, 2.9568, 5.4739, 5.4739, 5.4739, 0.]
    )

    def motion(q, v):
        return BodyMotion(mapping @ q, mapping, np.zeros(6))

    def damping_drag(_time, _q, v):
        body_speed = mapping @ v
        return -linear @ body_speed - drag * np.abs(body_speed) * body_speed

    dynamic = DynamicBody(
        rigid_mass=inertia,
        added_mass=(np.asarray(hydro["fAddedMass"]),),
        damping=(zero,),
        restoring=np.asarray(hydro["linearHydroRestCoef"]),
        static_force=np.array([
            0., 0., (rho * float(np.asarray(body.dispVol).item()) - 1156.5) * g,
            0., 0., 0.,
        ]),
        reference_position=np.r_[center, np.zeros(3)],
        motion=motion,
        excitation=lambda _time: np.zeros(6),
        radiation_kernel=np.asarray(hydro["irkb"]),
        state_excitation=damping_drag,
    )
    mooring_stiffness = np.diag([24_000., 5_000., 1_000.])
    mooring_damping = np.diag([5., 5., 0.])
    system = GeneralizedDynamics(
        (dynamic,), 3, pto_stiffness=mooring_stiffness,
        pto_damping=mooring_damping,
    )
    response = system.integrate(
        dt=dt, end_time=end_time, applied_force_history=applied,
    )
    q, v = response.coordinate, response.speed
    body_speed = response.body_velocity[:, 0, :]
    mooring_force = (-q @ mooring_stiffness.T
                     - v @ mooring_damping.T) @ mapping.T
    return WaveBotImpedanceResponse(
        time=response.time,
        body_position=response.body_position[:, 0, :],
        body_velocity=body_speed,
        body_acceleration=response.acceleration @ mapping.T,
        applied_force=applied,
        mooring_force=mooring_force,
        linear_damping_force=-body_speed @ linear.T,
        drag_force=-drag * np.abs(body_speed) * body_speed,
    )

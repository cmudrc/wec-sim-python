"""Published barge rigid and flexible trajectory against pinned MATLAB."""

import os
from pathlib import Path

import numpy as np
from scipy.io import loadmat

from wecsim import RegularWave, WEC


APPLICATIONS = os.environ.get("WEC_SIM_APPLICATIONS_DIR")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
LIVE = bool(APPLICATIONS and REFERENCE)
PORTABLE = Path(__file__).resolve().parent / "fixtures/matlab/gbm_barge.npz"
PORTABLE_HYDRO = Path(__file__).resolve().parent / "fixtures/matlab/gbm_barge.h5"


def _load(name):
    if LIVE:
        values = np.loadtxt(Path(REFERENCE) / f"GBM_BARGE_{name}.csv",
                            delimiter=",", ndmin=2)
        # Retain only the rigid time, pose, speed, and acceleration columns
        # used by both the live and compact-source checks.
        return values[:, np.r_[0:13, 43:49]] if name == "Barge_body1" else values
    with np.load(PORTABLE, allow_pickle=False) as source:
        return source[name]


def _error(actual, expected, limit, label):
    assert actual.shape == expected.shape, label
    assert np.isfinite(actual).all() and np.isfinite(expected).all(), label
    error = float(np.max(np.abs(actual - expected)))
    assert error < limit, f"{label}: {error:.6g} exceeds {limit}"


def test_published_barge_rigid_and_flexible_motion():
    if LIVE:
        flex = loadmat(Path(REFERENCE) / "GBM_BARGE_Flex_out.mat",
                       simplify_cells=True)["Flex_out"]
        time = np.asarray(flex["time"])
        values = np.asarray(flex["signals"]["values"])
        hydro = (Path(APPLICATIONS)
                 / "Generalized_Body_Modes/hydroData/barge.h5").resolve()
    else:
        with np.load(PORTABLE, allow_pickle=False) as source:
            time = source["time"]
            values = source["flex"]
        hydro = PORTABLE_HYDRO
    rigid = _load("Barge_body1")
    assert values.shape == (8001, 40)
    assert rigid.shape == (8001, 19)
    _error(time, rigid[:, 0], 1e-10, "flexible and rigid time")
    mode_q, mode_v, mode_a, total = np.split(values[:, :16], 4, axis=1)
    forces = values[:, 16:].reshape(-1, 6, 4)

    # The published Simulink block emits excitation, radiation, added mass,
    # restoring, viscous, and linear damping in that order.
    _error(total, forces[:, 0] - forces[:, 1:].sum(axis=1),
           1e-8, "source flexible force balance")
    mass = _load("flexible_effective_mass")
    stiffness = _load("flexible_stiffness")
    damping = _load("flexible_damping")
    _error(mode_a @ mass.T + mode_q @ stiffness.T + mode_v @ damping.T,
           total, 1e-6, "source flexible acceleration equation")
    q = np.column_stack((rigid[:, 1:7], mode_q))
    v = np.column_stack((rigid[:, 7:13], mode_v))
    restoring = _load("hydrostatic_stiffness")
    radiation = _load("radiation_damping")
    _error(forces[:, 1], v @ radiation[6:, :].T, 1e-6,
           "source flexible radiation")
    _error(forces[:, 3], q @ restoring[6:, :].T, 1e-6,
           "source flexible restoring")
    assert np.max(np.abs(forces[:, 4:])) == 0
    added = _load("added_mass")
    rigid_acceleration = rigid[:, 13:19]
    extrapolated = np.vstack((rigid_acceleration[:1], rigid_acceleration[:1],
                              2 * rigid_acceleration[1:-1]
                              - rigid_acceleration[:-2]))
    _error(forces[:, 2], extrapolated @ added[6:, :6].T, .05,
           "source sampled added-mass feedback")

    wec = WEC("GBM barge")
    barge = wec.body("barge", hydro, mass="equilibrium",
                     inertia=(6.667e7, 2.167e9, 2.167e9))
    wec.floating_gbm(barge)
    python = wec.run(RegularWave(height=2, period=8), dt=.05,
                     end_time=400, ramp_time=100)
    assert python.case["constraint"]["kind"] == "floating_gbm"
    modes = python.flexible_modes["barge"]
    _error(python.time, time, 1e-10, "paired time")
    _error(python.bodies["barge"].position, rigid[:, 1:7], 1.5e-3,
           "published rigid position")
    _error(python.bodies["barge"].velocity, rigid[:, 7:13], 1.2e-3,
           "published rigid velocity")
    _error(modes.position, mode_q, 6e-5,
           "published flexible displacement")
    _error(modes.velocity, mode_v, 5e-5,
           "published flexible velocity")
    _error(modes.acceleration, mode_a, 2e-4,
           "published flexible acceleration")
    _error(python.wave_elevation, _load("wave")[:, 1], 1e-12,
           "published regular wave")

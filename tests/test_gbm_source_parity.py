"""Pair generalized-body-mode preprocessing with the pinned MATLAB barge."""

import os
from pathlib import Path

import numpy as np

from wecsim.hydrodynamics.bodyClass import BodyClass


APPLICATIONS = os.environ.get("WEC_SIM_APPLICATIONS_DIR")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
FIXTURES = Path(__file__).parent / "fixtures/matlab"
FRESH = bool(APPLICATIONS and REFERENCE)


def _load(name):
    if FRESH:
        return np.loadtxt(Path(REFERENCE) / f"GBM_BARGE_{name}.csv",
                          delimiter=",", ndmin=2)
    archive = ("gbm_barge_preprocess.npz" if name == "excitation_coefficients"
               else "gbm_barge.npz")
    key = "coefficients" if name == "excitation_coefficients" else name
    with np.load(FIXTURES / archive) as saved:
        return saved[key]


def _logged_excitation():
    if FRESH:
        body = _load("Barge_body1")
        assert body.shape == (8001, 49)
        return body[:, 0], body[:, 19:25]
    with np.load(FIXTURES / "gbm_barge_preprocess.npz") as saved:
        return saved["body_time"], saved["rigid_excitation"]


def _max_error(actual, expected, limit, label):
    assert actual.shape == expected.shape, label
    assert np.isfinite(actual).all() and np.isfinite(expected).all(), label
    error = float(np.max(np.abs(actual - expected)))
    assert error < limit, f"{label}: {error:.6g} exceeds {limit}"


def test_barge_ten_dof_bem_and_flexible_matrices():
    hydro = ((Path(APPLICATIONS)
              / "Generalized_Body_Modes/hydroData/barge.h5").resolve()
             if FRESH else FIXTURES / "gbm_barge.h5")
    body = BodyClass(str(hydro))
    body.bodyNumber = body.bodyTotal = 1
    body.readH5file()
    assert int(body.dof.item()) == 10
    assert int(body.dof_gbm.item()) == 4
    body.mass = "equilibrium"
    body.hydroForcePre(
        2 * np.pi / 8, [0], 1, np.array([0.0]), [], .05, 1000, 9.81,
        "regular", np.zeros((2, 1)), 1, 1, 0, 0, 0,
    )
    force = body.hydroForce
    gbm = force["gbm"]
    _max_error(gbm["mass_ff"], _load("flexible_effective_mass"), 1e-6,
               "flexible effective mass")
    _max_error(gbm["stiffness"], _load("flexible_stiffness"), 1e-6,
               "flexible stiffness")
    _max_error(gbm["damping"], _load("flexible_damping"), 1e-9,
               "flexible mechanical damping")
    _max_error(force["linearHydroRestCoef"],
               _load("hydrostatic_stiffness"), 1e-5,
               "ten-DOF hydrostatic stiffness")
    _max_error(force["fAddedMass"], _load("added_mass"), 1e-5,
               "ten-DOF applied added mass")
    _max_error(force["fDamping"], _load("radiation_damping"), 1e-5,
               "ten-DOF radiation damping")
    coefficients = np.column_stack((force["fExt"]["re"],
                                    force["fExt"]["im"]))
    _max_error(coefficients, _load("excitation_coefficients"), 1e-5,
               "ten-DOF excitation coefficients")
    assert np.any(np.abs(coefficients[6:]) > 1), "flexible excitation missing"
    assert gbm["state_space"]["A"].shape == (8, 8)

    body_time, rigid_excitation = _logged_excitation()
    wave = _load("wave")
    assert body_time.shape == (8001,)
    assert rigid_excitation.shape == (8001, 6)
    assert wave.shape == (8001, 2)
    _max_error(body_time, wave[:, 0], 1e-10, "time grid")
    time = wave[:, 0]
    ramp = np.ones_like(time)
    early = time < 100
    ramp[early] = (1 - np.cos(np.pi * time[early] / 100)) / 2
    _max_error(wave[:, 1], ramp * np.cos(2 * np.pi * time / 8),
               1e-12, "regular wave")
    excitation = ramp[:, None] * (
        np.cos(2 * np.pi * time[:, None] / 8) * coefficients[:, 0]
        - np.sin(2 * np.pi * time[:, None] / 8) * coefficients[:, 1]
    )
    _max_error(excitation[:, :6], rigid_excitation, 1e-5,
               "logged rigid excitation")

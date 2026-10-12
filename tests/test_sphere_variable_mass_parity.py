"""Paired validation of the published variable-mass sphere application."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim import HydroState, RegularWave, WEC, WorldPoint
from wecsim.hydrodynamics.bodyClass import BodyClass


APPLICATIONS = os.environ.get("WEC_SIM_APPLICATIONS_DIR")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
pytestmark = pytest.mark.skipif(
    not (APPLICATIONS and REFERENCE),
    reason="paired MATLAB variable-mass output and Applications absent",
)


def _source():
    folder = Path(REFERENCE)
    body = np.loadtxt(folder / "SPHERE_VARIABLE_MASS_Variable_Mass_body1.csv",
                      delimiter=",")
    pto = np.loadtxt(folder / "SPHERE_VARIABLE_MASS_Variable_Mass_pto1.csv",
                     delimiter=",")
    state = np.loadtxt(folder / "SPHERE_VARIABLE_MASS_active_state.csv",
                       delimiter=",")
    properties = np.loadtxt(folder / "SPHERE_VARIABLE_MASS_states.csv",
                            delimiter=",")
    wave = np.loadtxt(folder / "SPHERE_VARIABLE_MASS_wave.csv",
                      delimiter=",")
    assert body.shape == (9001, 49)
    assert pto.shape == (9001, 25)
    assert state.shape == (9001, 2)
    assert properties.shape == (9, 6)
    assert wave.shape == (90001, 2)
    return body, pto, state, properties, wave


def _model(properties):
    hydro_dir = (Path(APPLICATIONS) / "Variable_Hydro" / "Variable_Mass"
                 / "hydroData")
    wec = WEC("Variable mass sphere")
    sphere = wec.variable_body(
        "sphere",
        [HydroState(hydro_dir / f"draft{i+1}.h5",
                    float(properties[i, 2]), tuple(properties[i, 3:6]))
         for i in range(9)],
        switch_times=range(100, 900, 100),
    )
    wec.coordinate("heave", sphere.move("heave"))
    wec.pto("PTO1", WorldPoint(0, 0, 2), sphere.at(0, 0, 0),
            axis=(0, 0, 1), damping=200_000)
    return wec


def _max_error(actual, expected, limit, label):
    assert actual.shape == expected.shape, label
    assert np.isfinite(actual).all() and np.isfinite(expected).all(), label
    error = float(np.max(np.abs(actual - expected)))
    assert error < limit, f"{label}: {error:.6g} exceeds {limit}"


def test_source_uses_all_nine_draft_states_and_matches_bem_forces():
    body, pto, active, properties, wave = _source()
    time = body[:, 0]
    expected_index = np.minimum(np.floor(time / 100).astype(int) + 1, 9)
    _max_error(active[:, 0], time, 1e-10, "active-state time")
    np.testing.assert_array_equal(active[:, 1], expected_index)
    _max_error(properties[:, 2],
               1025 * np.pi / 3 * np.arange(1, 10)**2
               * (15 - np.arange(1, 10)), 1e-6,
               "published displaced mass")
    _max_error(wave[:, 1], .5 * np.cos(2*np.pi*wave[:, 0]/8),
               1e-12, "source wave")

    coefficients = []
    for i in range(1, 10):
        path = (Path(APPLICATIONS) / "Variable_Hydro" / "Variable_Mass"
                / "hydroData" / f"draft{i}.h5")
        hydro = BodyClass(str(path))
        hydro.bodyNumber = hydro.bodyTotal = 1
        hydro.readH5file()
        hydro.mass = properties[i-1, 2]
        hydro.hydroStiffness = np.zeros((6, 6))
        hydro.viscDrag = {"Drag": np.zeros((6, 6)), "cd": np.zeros(6),
                          "characteristicArea": np.zeros(6)}
        hydro.linearDamping = np.zeros((6, 6))
        hydro.hydroForcePre(2*np.pi/8, [0], 1, np.array([0.0]), [],
                            .01, 1025, 9.81, "regular", np.zeros((2, 1)),
                            1, 1, 0, 0, 0)
        force = hydro.hydroForce
        coefficients.append((
            float(np.asarray(hydro.dispVol).item()),
            float(np.asarray(hydro.cg).ravel()[2]),
            force["fDamping"][2, 2],
            force["linearHydroRestCoef"][2, 2],
            force["fExt"]["re"][2], force["fExt"]["im"][2],
        ))
    values = np.asarray(coefficients)[expected_index - 1]
    volume, center, damping, stiffness, real, imaginary = values.T
    excitation = .5 * (real * np.cos(2*np.pi*time/8)
                       - imaginary * np.sin(2*np.pi*time/8))
    restoring = (stiffness * (body[:, 3] - center)
                 + (properties[expected_index-1, 2] - 1025*volume)*9.81)
    _max_error(excitation, body[:, 21], 1e-6, "state-selected excitation")
    _max_error(damping * body[:, 9], body[:, 27], 1e-6,
               "state-selected radiation damping")
    _max_error(restoring, body[:, 39], 1e-6,
               "state-selected restoring and weight")
    _max_error(body[:, 21] - body[:, 27] - body[:, 33]
               - body[:, 39], body[:, 15], 1e-6,
               "source hydrodynamic force assembly")
    _max_error(-200_000 * body[:, 9], pto[:, 15], 1e-6,
               "source PTO damping")


def test_public_python_variable_body_tracks_published_trajectory():
    body, pto, active, properties, wave = _source()
    result = _model(properties).run(
        RegularWave(height=1, period=8), dt=.01, end_time=900,
        ramp_time=0, rho=1025,
    )
    sample = slice(None, None, 10)
    extras = dict(result.raw.extra_outputs)
    _max_error(result.time[sample], body[:, 0], 1e-10, "time")
    _max_error(result.bodies["sphere"].position[sample, 2],
               body[:, 3], .003, "sphere heave")
    _max_error(result.bodies["sphere"].velocity[sample, 2],
               body[:, 9], .005, "sphere heave velocity")
    _max_error(result.ptos["PTO1"].stroke[sample], pto[:, 3],
               .003, "PTO stroke")
    _max_error(result.ptos["PTO1"].force[sample], pto[:, 15],
               1_000, "PTO force")
    _max_error(result.wave_elevation, wave[:, 1], 1e-12,
               "wave elevation")
    np.testing.assert_array_equal(extras["body1_hydro_state"][sample],
                                  active[:, 1])

"""Gate source force-clock attribution and an independent Python trajectory."""

import json
import os
from pathlib import Path

import numpy as np
import pytest

from wecsim import RegularCICWave, RegularWave, WEC, WorldPoint


ROOT = Path(__file__).resolve().parents[1]
RECORDS = os.environ.get("WEC_SIM_ELLIPSOID_CLOCK_DIR")
APPLICATIONS = os.environ.get("WEC_SIM_APPLICATIONS_DIR")
SELECTED_CASE = os.environ.get("WEC_SIM_ELLIPSOID_CASE")
FIXTURE = ROOT / "tests/fixtures/matlab/ellipsoid_ode45_nonlinear_clock/source_maxstep_common_grid.npz"
REFINED = ROOT / "tests/fixtures/matlab/ellipsoid_ode45_maxstep_only/refined_005_common_grid.npz"
pytestmark = pytest.mark.skipif(
    not (RECORDS and APPLICATIONS), reason="fresh pinned MATLAB clock audit absent",
)


def _read(folder, mode):
    shapes = {"body": (3001, 25), "pto": (3001, 19),
              "wave": (3001, 2), "settings": (1, 5)}
    record = {}
    for name, shape in shapes.items():
        values = np.loadtxt(folder / f"{mode}_{name}.csv", delimiter=",", ndmin=2)
        assert values.shape == shape and np.isfinite(values).all(), (mode, name)
        record[name] = values
    expected = [.005, .05, .005 if mode == "fine" else .05,
                .05 if mode == "coarse_force_cic" else .005, .005]
    np.testing.assert_allclose(record["settings"][0], expected, rtol=0, atol=1e-12)
    time = np.arange(3001) * .05
    for name in ("body", "pto", "wave"):
        np.testing.assert_allclose(record[name][:, 0], time, rtol=0, atol=1e-10)
    return record


def _error(a, b):
    return float(np.max(np.abs(a - b)))


@pytest.mark.parametrize(
    "wave_case", [SELECTED_CASE] if SELECTED_CASE else ["Regular", "RegularCIC"],
)
def test_nonlinear_and_cic_clocks_separately_explain_maxstep_gap(wave_case):
    key = "regular" if wave_case == "Regular" else "regular_cic"
    folder = Path(RECORDS) / wave_case
    fine = _read(folder, "fine")
    coarse = _read(folder, "coarse_force")
    both_coarse = (_read(folder, "coarse_force_cic")
                   if wave_case == "RegularCIC" else coarse)
    with (np.load(FIXTURE, allow_pickle=False) as prior,
          np.load(REFINED, allow_pickle=False) as refined):
        source = {name: prior[f"{key}_maxstep005_{name}"]
                  for name in ("time", "heave", "speed", "restoring",
                               "excitation", "pto_force", "wave")}
        fine_reference = {name: refined[f"{key}_{name}"]
                          for name in ("heave", "speed", "restoring",
                                       "pto_force", "wave")}

    # All arms inherit the same wave. The regenerated fine-clock source also
    # matches the archived fully refined source at every common sample.
    for record in (fine, coarse, both_coarse):
        assert _error(record["wave"][:, 1], source["wave"]) < 1e-10
    assert _error(fine["body"][:, 3], fine_reference["heave"]) < 1e-7
    assert _error(fine["body"][:, 9], fine_reference["speed"]) < 1e-7
    assert _error(fine["body"][:, 15], fine_reference["restoring"]) < .1
    assert _error(fine["pto"][:, 15], fine_reference["pto_force"]) < .1

    force_clock_heave = _error(fine["body"][:, 3], coarse["body"][:, 3])
    force_clock_speed = _error(fine["body"][:, 9], coarse["body"][:, 9])
    force_clock_restoring = _error(fine["body"][:, 15], coarse["body"][:, 15])
    assert .015 < force_clock_heave < .03
    assert .02 < force_clock_speed < .04
    assert 30_000 < force_clock_restoring < 60_000

    # Restoring both source force clocks reproduces the MaxStep-only arm,
    # which kept the published 0.05 s force clocks. This comparison isolates
    # sample timing; it makes no claim that the published solve is stable.
    errors = {
        "heave": _error(both_coarse["body"][:, 3], source["heave"]),
        "speed": _error(both_coarse["body"][:, 9], source["speed"]),
        "restoring": _error(both_coarse["body"][:, 15], source["restoring"]),
        "excitation": _error(both_coarse["body"][:, 21], source["excitation"]),
        "pto_force": _error(both_coarse["pto"][:, 15], source["pto_force"]),
    }
    assert errors["heave"] < 1e-7 and errors["speed"] < 1e-7
    assert max(errors["restoring"], errors["excitation"],
               errors["pto_force"]) < .1
    if wave_case == "RegularCIC":
        cic_heave = _error(coarse["body"][:, 3], both_coarse["body"][:, 3])
        cic_speed = _error(coarse["body"][:, 9], both_coarse["body"][:, 9])
        assert 5e-5 < cic_heave < 2e-4
        assert 2e-4 < cic_speed < 8e-4
    else:
        cic_heave = cic_speed = 0.0

    # The public physical solve uses only pinned HDF5 and wave settings. Its
    # proximity to the fine-clock source links this attribution to parity.
    app = Path(APPLICATIONS) / "Nonlinear_Hydro"
    wec = WEC("ellipsoid source-clock audit")
    body = wec.body(
        "ellipsoid", app / "hydroData/ellipsoid.h5", mass="equilibrium",
        inertia=(1.375264e6, 1.375264e6, 1.341721e6),
        geometry_file=app / "geometry/elipsoid.stl",
        nonlinear_hydro="instantaneous", drag_coefficient=1,
        drag_area=np.pi * 25,
    )
    wec.coordinate("heave", body.move("heave"))
    wec.pto("PTO1", WorldPoint(0, 0, -12.5), body.at(0, 0, 0),
            damping=1_200_000)
    wave = RegularWave(4, 6) if wave_case == "Regular" else RegularCICWave(4, 6)
    python = wec.run(
        wave, dt=.05, end_time=150, ramp_time=50, rho=1025,
        radiation_memory=60 if wave_case == "RegularCIC" else None,
    )
    python_heave = python.bodies["ellipsoid"].position[:, 2]
    python_to_fine = _error(python_heave, fine["body"][:, 3])
    python_to_coarse = _error(python_heave, coarse["body"][:, 3])
    assert python_to_fine < (.003 if wave_case == "Regular" else .0035)
    assert python_to_coarse > 3 * python_to_fine
    print("ellipsoid clock attribution:", json.dumps({
        "case": wave_case, "force_clock_heave_m": force_clock_heave,
        "force_clock_speed_mps": force_clock_speed,
        "force_clock_restoring_n": force_clock_restoring,
        "cic_clock_heave_m": cic_heave, "cic_clock_speed_mps": cic_speed,
        "restored_to_maxstep": errors,
        "python_to_fine_heave_m": python_to_fine,
        "python_to_coarse_heave_m": python_to_coarse,
    }, sort_keys=True))

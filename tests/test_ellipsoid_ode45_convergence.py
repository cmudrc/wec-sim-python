"""Full-duration independent motion against a source-step-screened ode45 input.

The unchanged 0.05 s published cases retain their separate diagnostic gates.
This test concerns copied, derived source inputs at 0.005, 0.0025, 0.00125,
and 0.000625 s.
"""

import os
from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import trapezoid

from wecsim import RegularCICWave, RegularWave, WEC, WorldPoint
from wecsim.hydrodynamics.nonlinearHydro import HeaveMeshHydro


ROOT = Path(__file__).resolve().parents[1]
APPLICATIONS = os.environ.get("WEC_SIM_APPLICATIONS_DIR")
RECORDS = os.environ.get("WEC_SIM_ELLIPSOID_CONVERGENCE_DIR")
CASE = os.environ.get("WEC_SIM_ELLIPSOID_CASE")
ARCHIVED_005 = (ROOT / "tests/fixtures/matlab/ellipsoid_ode45_maxstep_only"
                / "refined_005_common_grid.npz")
pytestmark = pytest.mark.skipif(
    not (APPLICATIONS and RECORDS and CASE),
    reason="fresh pinned R2025b ode45 convergence output absent",
)

# These original 150 s physical limits come from the established refined-step
# gate. They are deliberately unchanged; the source screen adds an independent
# requirement that the final source change be below one third of these limits.
PYTHON_LIMITS = {
    "Regular": (.003, .004, 5000, 4000),
    "RegularCIC": (.0035, .0045, 6000, 5000),
}
SOURCE_STABILITY_LIMITS = {
    "Regular": (.001, .0013, 1600, 1300),
    "RegularCIC": (.0011, .0015, 2000, 1600),
}
PYTHON_STEP_LIMITS = (.0001, .0001, 100, 100)
# The source signed-work sequence must settle independently before its finest
# member can serve as the derived reference for the Python 0.1% work gate.
SOURCE_WORK_STEP_LIMIT = .001
STEPS = {"step005": .005, "step0025": .0025,
         "step00125": .00125, "step000625": .000625}
# Columns: time, heave, speed, PTO stroke, stroke speed, PTO force,
# absorbed power, wave, restoring, excitation, acceleration.
CHANNELS = ((1, "heave position"), (2, "heave speed"),
            (5, "PTO force"), (6, "absorbed power"))


def _error(a, b):
    assert a.shape == b.shape
    assert np.isfinite(a).all() and np.isfinite(b).all()
    return float(np.max(np.abs(a - b)))


def _check(a, b, limit, label):
    error = _error(a, b)
    print(f"{label}: {error:.8g} (gate {limit:.8g})")
    assert error < limit, f"{label}: {error:.8g} >= {limit:.8g}"
    return error


def _read_source(case):
    folder = Path(RECORDS) / case
    result = {}
    for label, step in STEPS.items():
        history = np.loadtxt(folder / f"{label}_history.csv",
                             delimiter=",", ndmin=2)
        settings = np.loadtxt(folder / f"{label}_settings.csv",
                              delimiter=",", ndmin=2)
        assert history.shape == (round(150 / step) + 1, 11), label
        assert np.isfinite(history).all(), label
        np.testing.assert_allclose(settings, np.full((1, 5), step),
                                   rtol=0, atol=1e-12)
        np.testing.assert_allclose(history[:, 0],
                                   np.arange(len(history)) * step,
                                   rtol=0, atol=1e-9)
        result[label] = history
    return result


def _python(case, step, app):
    wec = WEC("source-step-screened ellipsoid ode45")
    body = wec.body(
        "ellipsoid", app / "hydroData/ellipsoid.h5",
        mass="equilibrium", inertia=(1.375264e6, 1.375264e6, 1.341721e6),
        geometry_file=app / "geometry/elipsoid.stl",
        nonlinear_hydro="instantaneous", drag_coefficient=1,
        drag_area=np.pi * 25,
    )
    wec.coordinate("heave", body.move("heave"))
    wec.pto("PTO1", WorldPoint(0, 0, -12.5), body.at(0, 0, 0),
            damping=1_200_000)
    wave = RegularCICWave(4, 6) if case == "RegularCIC" else RegularWave(4, 6)
    return wec.run(wave, dt=step, end_time=150, ramp_time=50,
                   radiation_memory=60 if case == "RegularCIC" else None,
                   rho=1025)


def _python_channels(result):
    body = result.bodies["ellipsoid"]
    pto = result.ptos["PTO1"]
    return (body.position[:, 2], body.velocity[:, 2],
            pto.force, pto.absorbed_power)


def _validate_full_gate(source, case, app):
    assert case in PYTHON_LIMITS
    for label, step in STEPS.items():
        history = source[label]
        assert history.shape == (round(150 / step) + 1, 11)
        assert np.isfinite(history).all()
        np.testing.assert_allclose(history[:, 0],
                                   np.arange(len(history)) * step,
                                   rtol=0, atol=1e-9)
    five, two_half, one_quarter, finest = (
        source[label] for label in STEPS
    )
    key = "regular_cic" if case == "RegularCIC" else "regular"
    with np.load(ARCHIVED_005, allow_pickle=False) as archived:
        for column, suffix, limit in ((0, "time", 1e-10),
                                      (1, "heave", 1e-7),
                                      (2, "speed", 1e-7),
                                      (5, "pto_force", .1),
                                      (7, "wave", 1e-10),
                                      (8, "restoring", .1)):
            _check(five[::10, column], archived[f"{key}_{suffix}"],
                   limit, f"regenerated 0.005 s {suffix}")

    _check(five[:, 7], two_half[::2, 7], 1e-10,
           "source 0.005/0.0025 s wave identity")
    _check(two_half[:, 7], one_quarter[::2, 7], 1e-10,
           "source 0.0025/0.00125 s wave identity")
    _check(one_quarter[:, 7], finest[::2, 7], 1e-10,
           "source 0.00125/0.000625 s wave identity")
    for (column, channel), stability in zip(
        CHANNELS, SOURCE_STABILITY_LIMITS[case], strict=True
    ):
        first = _error(five[:, column], two_half[::2, column])
        previous = _error(two_half[:, column], one_quarter[::2, column])
        final = _check(one_quarter[:, column], finest[::2, column],
                       stability, f"source final-step {channel}")
        assert previous < .7 * first, (
            f"source {channel} did not contract: {first:.8g} -> {previous:.8g}"
        )
        assert final < .7 * previous, (
            f"source {channel} did not contract: {previous:.8g} -> {final:.8g}"
        )

    source_works = [float(trapezoid(history[:, 6], history[:, 0]))
                    for history in (five, two_half, one_quarter, finest)]
    assert all(np.isfinite(work) and work > 1e3 for work in source_works)
    work_changes = [abs(later - earlier) / later
                    for earlier, later in zip(source_works, source_works[1:])]
    print("source signed PTO work (J):", source_works)
    print("source signed PTO work relative step changes:", work_changes)
    assert work_changes[-1] < SOURCE_WORK_STEP_LIMIT
    assert work_changes[1] < .7 * work_changes[0]
    assert work_changes[2] < .7 * work_changes[1]

    mesh = HeaveMeshHydro.from_stl(
        app / "geometry/elipsoid.stl", center_z=-2, rho=1025,
        gravity=9.81, depth=70, period=6, height=4, ramp_time=50,
        mass=None, drag_coefficient=1, drag_area=np.pi * 25,
    )
    samples = np.linspace(1, len(finest) - 1, 1001, dtype=int)
    previous = np.array([
        mesh.forces(finest[i - 1, 0], finest[i - 1, 1] + 2,
                    finest[i - 1, 2])[0]
        for i in samples
    ])
    current = np.array([
        mesh.forces(finest[i, 0], finest[i, 1] + 2,
                    finest[i, 2])[0]
        for i in samples
    ])
    _check(previous, -finest[samples, 8], 1e-6,
           "finest source preceding-state restoring law")
    assert _error(current, -finest[samples, 8]) > 500

    python_two_half = _python(case, .0025, app)
    python_one_quarter = _python(case, .00125, app)
    _check(python_two_half.time, python_one_quarter.time[::2], 1e-9,
           "Python common times")
    for (_, name), coarse, fine, limit in zip(
        CHANNELS, _python_channels(python_two_half),
        _python_channels(python_one_quarter), PYTHON_STEP_LIMITS, strict=True
    ):
        _check(coarse, fine[::2], limit, f"Python final-step {name}")
    _check(python_one_quarter.time, finest[::2, 0], 1e-9,
           "Python/source time")
    _check(python_one_quarter.wave_elevation, finest[::2, 7], 1e-10,
           "Python/source wave")
    for (column, name), actual, limit in zip(
        CHANNELS, _python_channels(python_one_quarter),
        PYTHON_LIMITS[case], strict=True
    ):
        _check(actual, finest[::2, column], limit,
               f"independent physical {name}")
    _check(python_one_quarter.ptos["PTO1"].stroke,
           finest[::2, 3], PYTHON_LIMITS[case][0],
           "independent physical PTO stroke")
    source_work = source_works[-1]
    python_work = float(trapezoid(
        python_one_quarter.ptos["PTO1"].absorbed_power,
        python_one_quarter.time,
    ))
    assert source_work > 1e3
    work_error = abs(python_work / source_work - 1)
    print(f"independent signed PTO work: {work_error:.6%} (gate 0.1%)")
    assert work_error < .001


def test_pinned_source_step_screen_and_full_duration_independent_motion():
    assert CASE in PYTHON_LIMITS
    _validate_full_gate(_read_source(CASE), CASE,
                        Path(APPLICATIONS) / "Nonlinear_Hydro")

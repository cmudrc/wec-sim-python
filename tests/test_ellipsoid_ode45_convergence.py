"""Full-duration independent motion against a source-step-screened ode45 input.

The unchanged 0.05 s published cases retain their separate diagnostic gates.
This test concerns copied, derived source inputs at 0.005/0.0025/0.00125 s.
"""

import os
from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import trapezoid

from wecsim import RegularCICWave, RegularWave, WEC, WorldPoint


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
STEPS = {"step005": .005, "step0025": .0025,
         "step00125": .00125}
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


def test_pinned_source_step_screen_and_full_duration_independent_motion():
    assert CASE in PYTHON_LIMITS
    source = _read_source(CASE)
    five, two_half, one_quarter = (
        source[label] for label in STEPS
    )
    key = "regular_cic" if CASE == "RegularCIC" else "regular"
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
    for (column, channel), stability in zip(
        CHANNELS, SOURCE_STABILITY_LIMITS[CASE], strict=True
    ):
        previous = _error(five[:, column], two_half[::2, column])
        final = _check(two_half[:, column], one_quarter[::2, column],
                       stability, f"source final-step {channel}")
        assert final < .7 * previous, (
            f"source {channel} did not contract: {previous:.8g} -> {final:.8g}"
        )

    app = Path(APPLICATIONS) / "Nonlinear_Hydro"
    python_two_half = _python(CASE, .0025, app)
    python_one_quarter = _python(CASE, .00125, app)
    _check(python_two_half.time, python_one_quarter.time[::2], 1e-9,
           "Python common times")
    for (_, name), coarse, fine, limit in zip(
        CHANNELS, _python_channels(python_two_half),
        _python_channels(python_one_quarter), PYTHON_STEP_LIMITS, strict=True
    ):
        _check(coarse, fine[::2], limit, f"Python final-step {name}")
    _check(python_one_quarter.time, one_quarter[:, 0], 1e-9,
           "Python/source time")
    _check(python_one_quarter.wave_elevation, one_quarter[:, 7], 1e-10,
           "Python/source wave")
    for (column, name), actual, limit in zip(
        CHANNELS, _python_channels(python_one_quarter),
        PYTHON_LIMITS[CASE], strict=True
    ):
        _check(actual, one_quarter[:, column], limit,
               f"independent physical {name}")
    _check(python_one_quarter.ptos["PTO1"].stroke,
           one_quarter[:, 3], PYTHON_LIMITS[CASE][0],
           "independent physical PTO stroke")
    source_work = float(trapezoid(one_quarter[:, 6], one_quarter[:, 0]))
    python_work = float(trapezoid(
        python_one_quarter.ptos["PTO1"].absorbed_power,
        python_one_quarter.time,
    ))
    assert source_work > 1e3
    work_error = abs(python_work / source_work - 1)
    print(f"independent signed PTO work: {work_error:.6%} (gate 0.1%)")
    assert work_error < .001

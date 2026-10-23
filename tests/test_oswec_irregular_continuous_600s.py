"""Gate the derived continuous-heading OSWEC control against direct MATLAB."""

import json
from pathlib import Path

import numpy as np

from wecsim import PMWave, WEC, WorldPoint


FIXTURES = Path(__file__).parent / "fixtures/matlab"
CASE = FIXTURES / "oswec_irregular_continuous_600s"
HYDRO = FIXTURES / "oswec_multi_wave/oswec.h5"


def _max_gap(actual, expected, limit, label):
    gap = float(np.max(np.abs(actual - expected)))
    assert gap < limit, f"{label}: {gap} >= {limit}"


def _relative_work_gap(a, b):
    return abs(a - b) / abs(b)


def test_direct_matlab_continuous_source_is_step_stable():
    """The modified source model must pass before it can be a baseline."""
    limits = json.loads((CASE / "gate_spec.json").read_text())
    with np.load(CASE / "source.npz") as saved:
        coarse, fine = saved["coarse"], saved["fine"]
        coarse_wave, fine_wave = saved["coarse_wave"], saved["fine_wave"]
    assert coarse.shape == (60_001, 5)
    assert fine.shape == (120_001, 5)
    assert coarse_wave.shape == (60_001, 2)
    assert fine_wave.shape == (120_001, 2)
    assert np.isfinite(coarse).all() and np.isfinite(fine).all()
    for trace in (coarse, fine):
        _max_gap(trace[:, 3], -120_000 * trace[:, 2], 1e-5,
                 "source PTO damping law")
        _max_gap(trace[:, 4], trace[:, 3] * trace[:, 2], 1e-8,
                 "source-signed PTO power")
    _max_gap(coarse_wave[:, 0], coarse[:, 0], 1e-10, "coarse wave times")
    _max_gap(fine_wave[:, 0], fine[:, 0], 1e-10, "fine wave times")
    _max_gap(fine[::2, 0], coarse[:, 0], 1e-10, "source sample times")
    _max_gap(fine_wave[::2, 1], coarse_wave[:, 1],
             limits["wave_elevation_max_abs_m"], "source wave")
    for column, key, label in (
        (1, "yaw_max_abs_rad", "source yaw step change"),
        (2, "yaw_speed_max_abs_rad_s", "source speed step change"),
        (3, "pto_torque_max_abs_Nm", "source torque step change"),
    ):
        _max_gap(fine[::2, column], coarse[:, column], limits[key], label)
    coarse_work = -np.trapezoid(coarse[:, 4], coarse[:, 0])
    fine_work = -np.trapezoid(fine[:, 4], fine[:, 0])
    assert coarse_work > 0 and fine_work > 0
    assert (_relative_work_gap(fine_work, coarse_work)
            < limits["relative_pto_work_max"])


def test_physical_public_wec_matches_full_direct_matlab_continuous_case():
    """Advance the physical model from only the HDF5 and PM wave settings."""
    limits = json.loads((CASE / "gate_spec.json").read_text())
    with np.load(CASE / "source.npz") as saved:
        source, source_wave = saved["fine"], saved["fine_wave"]

    wec = WEC("OSWEC 600 s continuous-heading control")
    flap = wec.body("flap", HYDRO, mass=12_700, inertia=(1.85e6,) * 3,
                    passive_yaw=True, passive_yaw_threshold=0)
    wec.body("base", HYDRO, mass=999, inertia=(999,) * 3)
    yaw = wec.coordinate("yaw", flap.move("yaw", pivot=WorldPoint(0, 0, -8.9)))
    wec.rotational_pto("hinge", yaw, damping=120_000)
    result = wec.run(
        PMWave(2.5, 8, direction=10, seed=1, phase_generator="matlab"),
        dt=0.005, end_time=600, ramp_time=100, radiation_memory=40,
    )
    _max_gap(result.time, source[:, 0], 1e-10, "time")
    _max_gap(result.wave_elevation, source_wave[:, 1],
             limits["wave_elevation_max_abs_m"], "incident wave")
    angle = result.bodies["flap"].position[:, 5]
    speed = result.bodies["flap"].velocity[:, 5]
    torque = result.ptos["hinge"].force
    for actual, expected, key, label in (
        (angle, source[:, 1], "yaw_max_abs_rad", "flap yaw"),
        (speed, source[:, 2], "yaw_speed_max_abs_rad_s", "flap speed"),
        (torque, source[:, 3], "pto_torque_max_abs_Nm", "PTO torque"),
    ):
        _max_gap(actual, expected, limits[key], label)
    source_work = -np.trapezoid(source[:, 4], source[:, 0])
    python_work = np.trapezoid(result.ptos["hinge"].absorbed_power, result.time)
    assert source_work > 0 and python_work > 0
    assert (_relative_work_gap(python_work, source_work)
            < limits["relative_pto_work_max"])

    # This case evaluates a continuous heading; no discrete bank events are
    # prescribed or replayed. Both independent motions stay inside the BEM span.
    heading_bounds = limits["relative_heading_range_deg"]
    for yaw_trace in (angle, source[:, 1]):
        relative_heading = 10 - np.degrees(yaw_trace)
        assert relative_heading.min() > heading_bounds[0]
        assert relative_heading.max() < heading_bounds[1]

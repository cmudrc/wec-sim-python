"""Pair a second OSWEC regular sea through the public source-step solver."""

import hashlib
import json
from pathlib import Path

import numpy as np

from wecsim import RegularWave, WEC, WorldPoint


ROOT = Path(__file__).resolve().parents[1]
HYDRO = ROOT / "examples/data/oswec.h5"
FIXTURE = ROOT / "tests/fixtures/matlab/oswec_regular_derived_60s"
HYDRO_SHA256 = "28aa028d404a88f877924bc435e7ad3a4c50e26f9706c1f39ef7324e5ec60750"


def _gate(actual, expected, limit, label):
    assert actual.shape == expected.shape, label
    assert np.isfinite(actual).all() and np.isfinite(expected).all(), label
    error = float(np.max(np.abs(actual - expected)))
    assert error < limit, f"{label}: {error:.9g} exceeds {limit:.9g}"


def test_derived_regular_sea_full_motion_pto_and_work():
    """Change wave height and period, then advance independently for 60 s."""
    assert hashlib.sha256(HYDRO.read_bytes()).hexdigest() == HYDRO_SHA256
    spec = json.loads((FIXTURE / "gate_spec.json").read_text())
    params = spec["derived_input"]
    limits = spec["limits"]
    wec = WEC("OSWEC derived regular sea")
    flap = wec.body("flap", HYDRO, mass=127_000, inertia=(0, 1.85e6, 0))
    base = wec.body("base", HYDRO, mass=999, inertia=(999, 999, 999))
    wec.fixed_hinge(
        flap, base, location=WorldPoint(0, 0, -10),
        pto_location=WorldPoint(0, 0, -8.9), damping=12_000,
        added_mass_scheme="simulink_delay",
    )
    result = wec.run(
        RegularWave(params["wave_height_m"], params["wave_period_s"]),
        dt=params["step_s"], end_time=params["duration_s"],
        ramp_time=params["ramp_s"], solver="wecsim_ode4",
    )

    with np.load(FIXTURE / "source.npz", allow_pickle=False) as source:
        body, pto, wave = source["flap"], source["pto"], source["wave"]
    assert body.shape == (601, 7)
    assert pto.shape == (601, 5)
    assert wave.shape == (601, 2)
    np.testing.assert_allclose(result.time, body[:, 0], rtol=0, atol=1e-10)
    np.testing.assert_allclose(result.time, pto[:, 0], rtol=0, atol=1e-10)
    np.testing.assert_allclose(result.time, wave[:, 0], rtol=0, atol=1e-10)
    np.testing.assert_allclose(pto[:, 3], -12_000 * pto[:, 2],
                               rtol=0, atol=1e-5)

    flap_history = result.bodies["flap"]
    pto_history = result.ptos["hinge"]
    for key, actual, expected in (
        ("wave_elevation_m", result.wave_elevation, wave[:, 1]),
        ("flap_surge_m", flap_history.position[:, 0], body[:, 1]),
        ("flap_heave_m", flap_history.position[:, 2], body[:, 2]),
        ("flap_pitch_rad", flap_history.position[:, 4], body[:, 3]),
        ("flap_pitch_speed_rad_s", flap_history.velocity[:, 4], body[:, 6]),
        ("pto_stroke_rad", pto_history.stroke, pto[:, 1]),
        ("pto_speed_rad_s", pto_history.velocity, pto[:, 2]),
        ("pto_torque_Nm", pto_history.force, pto[:, 3]),
    ):
        _gate(actual, expected, limits[key], key)

    for start, stop in spec["work_windows_s"]:
        window = (result.time >= start) & (result.time <= stop)
        source_work = -np.trapezoid(pto[window, 4], result.time[window])
        python_work = np.trapezoid(
            pto_history.absorbed_power[window], result.time[window],
        )
        assert source_work > 0
        relative_error = abs(python_work - source_work) / source_work
        assert relative_error < limits["pto_work_relative"], (
            start, stop, relative_error,
        )

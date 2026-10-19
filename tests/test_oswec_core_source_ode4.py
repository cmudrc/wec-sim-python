"""Gate the published OSWEC PTO work through the opt-in source solver."""

from pathlib import Path

import numpy as np
import pytest

from wecsim import PMWave, WEC, WorldPoint


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "tests/fixtures/matlab/oswec_irregular_wave.npz"
HYDRO = ROOT / "examples/data/oswec.h5"
WINDOWS = ((0, 400), (0, 100), (100, 250), (250, 400))


def _work(time, power, start, stop):
    selected = (time >= start - 1e-9) & (time <= stop + 1e-9)
    return float(np.trapezoid(power[selected], time[selected]))


def _wec(added_mass_scheme="simulink_delay"):
    wec = WEC("Published OSWEC source ode4")
    flap = wec.body("flap", HYDRO, mass=127_000, inertia=(1.85e6,) * 3)
    base = wec.body("base", HYDRO, mass=999, inertia=(999,) * 3)
    wec.fixed_hinge(
        flap, base, location=WorldPoint(0, 0, -10),
        pto_location=WorldPoint(0, 0, -8.9), damping=12_000,
        added_mass_scheme=added_mass_scheme,
    )
    return wec


def _wave(source, phase_file):
    np.savetxt(phase_file, source["wave_components"][:, 3:], delimiter=",")
    directions = source["wave_directions"]
    return PMWave(
        2.5, 8, directions=tuple(directions[:, 0]),
        spreading=tuple(directions[:, 1]), phase_file=phase_file,
    )


def test_published_oswec_motion_and_pto_work_on_source_ode4(tmp_path):
    """Run the full 400 s sea without prescribed MATLAB motion or force."""
    with np.load(SOURCE, allow_pickle=False) as source:
        expected = {name: source[name].copy() for name in (
            "time", "wave_elevation", "body_position", "body_velocity",
            "pto_stroke", "pto_speed", "pto_force", "pto_power",
        )}
        wave = _wave(source, tmp_path / "phases.csv")

    result = _wec().run(
        wave, dt=.1, end_time=400, ramp_time=100,
        radiation_memory=30, solver="wecsim_ode4",
    )
    np.testing.assert_allclose(result.time, expected["time"], rtol=0,
                               atol=1e-10)
    np.testing.assert_allclose(result.wave_elevation,
                               expected["wave_elevation"], rtol=0,
                               atol=1e-10)
    np.testing.assert_allclose(result.bodies["flap"].position,
                               expected["body_position"][0], rtol=0,
                               atol=1e-8)
    np.testing.assert_allclose(result.bodies["flap"].velocity,
                               expected["body_velocity"][0], rtol=0,
                               atol=1e-8)
    pto = result.ptos["hinge"]
    for actual, name, bound in (
        (pto.stroke, "pto_stroke", 1e-8),
        (pto.velocity, "pto_speed", 1e-8),
        (pto.force, "pto_force", 1e-5),
        (pto.absorbed_power, "pto_power", 1e-6),
    ):
        expected_trace = expected[name]
        if name == "pto_power":
            expected_trace = -expected_trace
        np.testing.assert_allclose(actual, expected_trace, rtol=0, atol=bound)
    for start, stop in WINDOWS:
        source_work = _work(expected["time"], -expected["pto_power"],
                            start, stop)
        python_work = _work(result.time, pto.absorbed_power, start, stop)
        assert source_work > 1_000
        assert abs(python_work - source_work) / source_work < 1e-8


def test_source_ode4_requires_simulink_added_mass(tmp_path):
    with np.load(SOURCE, allow_pickle=False) as source:
        wave = _wave(source, tmp_path / "phases.csv")
    with pytest.raises(ValueError, match="simulink_delay"):
        _wec("implicit").run(
            wave, dt=.1, end_time=.1, ramp_time=100,
            radiation_memory=30, solver="wecsim_ode4",
        )

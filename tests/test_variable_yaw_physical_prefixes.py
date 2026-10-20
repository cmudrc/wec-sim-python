"""Gate OSWEC physical motion prefixes with the source-stability caveat."""

import json
from pathlib import Path

import numpy as np

from wecsim import PMWave, WEC, WorldPoint


FIXTURES = Path(__file__).parent / "fixtures/matlab"
# This pinned common OSWEC BEMIO export has the required 40 s radiation IRF.
HYDRO = FIXTURES / "oswec_multi_wave/oswec.h5"


def _selected_heading_wec(label):
    headings = np.arange(-30, 30.25, .25)
    wec = WEC(label)
    flap = wec.body("flap", HYDRO, mass=12_700, inertia=(1.85e6,) * 3,
                    passive_yaw=True, yaw_heading_bank=headings)
    wec.body("base", HYDRO, mass=999, inertia=(999,) * 3)
    yaw = wec.coordinate("yaw", flap.move("yaw", pivot=WorldPoint(0, 0, -8.9)))
    wec.rotational_pto("hinge", yaw, damping=120_000)
    return wec


def test_physical_python_matches_stable_140s_selected_heading_prefix():
    """Independently gate the source-stable prefix, not the 600 s input."""
    with np.load(FIXTURES / "oswec_yaw_160_source.npz") as saved:
        source = saved["quarter"][:56_001]
        source_wave = saved["quarter_wave"][:56_001]
    gates = json.loads((FIXTURES / "oswec_yaw_160_gate_spec.json").read_text())
    limits = gates["independent_python_source_limits"]

    wec = _selected_heading_wec("OSWEC selected-heading stable prefix")
    result = wec.run(
        PMWave(2.5, 8, direction=10, seed=1, phase_generator="matlab"),
        dt=.0025, end_time=140, ramp_time=100, radiation_memory=40,
    )
    np.testing.assert_allclose(result.time, source[:, 0], rtol=0, atol=1e-10)
    np.testing.assert_allclose(result.wave_elevation, source_wave[:, 1],
                               rtol=0, atol=1e-12)

    source_work = -np.trapezoid(source[:, 6], source[:, 0])
    python_work = np.trapezoid(result.ptos["hinge"].absorbed_power, result.time)
    assert source_work > 0 and python_work > 0
    measured = {
        "yaw_rad": np.max(np.abs(result.bodies["flap"].position[:, 5] - source[:, 1])),
        "yaw_speed_rad_s": np.max(np.abs(result.bodies["flap"].velocity[:, 5] - source[:, 2])),
        "pto_torque_Nm": np.max(np.abs(result.ptos["hinge"].force - source[:, 5])),
        "relative_pto_work": abs(python_work - source_work) / source_work,
    }
    for channel, error in measured.items():
        assert error < limits[channel], f"{channel}: {error:.6g} >= {limits[channel]}"


def test_finer_physical_python_matches_160s_quarter_step_diagnostic():
    """Extend the independent comparison without promoting the unstable source case.

    MATLAB's half-to-quarter-step self-change fails these same limits at 160 s,
    so passing this one quarter-step trace is only a bounded diagnostic.
    """
    with np.load(FIXTURES / "oswec_yaw_160_source.npz", allow_pickle=False) as saved:
        source = saved["quarter"]
        source_wave = saved["quarter_wave"]
    gates = json.loads((FIXTURES / "oswec_yaw_160_gate_spec.json").read_text())
    assert gates["published_case_promotion"] is False
    limits = gates["independent_python_source_limits"]

    wec = _selected_heading_wec("OSWEC selected-heading finer physical diagnostic")
    result = wec.run(
        PMWave(2.5, 8, direction=10, seed=1, phase_generator="matlab"),
        dt=.00125, end_time=160, ramp_time=100, radiation_memory=40,
    )

    # The source trace is sampled at .0025 s. Python advances twice per
    # source sample and has no source body state, force, or heading-event input.
    sampled = slice(None, None, 2)
    np.testing.assert_allclose(result.time[sampled], source[:, 0],
                               rtol=0, atol=1e-10)
    np.testing.assert_allclose(result.wave_elevation[sampled], source_wave[:, 1],
                               rtol=0, atol=1e-12)
    source_work = -np.trapezoid(source[:, 6], source[:, 0])
    python_work = np.trapezoid(result.ptos["hinge"].absorbed_power, result.time)
    assert source_work > 0 and python_work > 0
    measured = {
        "yaw_rad": np.max(np.abs(result.bodies["flap"].position[sampled, 5] - source[:, 1])),
        "yaw_speed_rad_s": np.max(np.abs(result.bodies["flap"].velocity[sampled, 5] - source[:, 2])),
        "pto_torque_Nm": np.max(np.abs(result.ptos["hinge"].force[sampled] - source[:, 5])),
        "relative_pto_work": abs(python_work - source_work) / source_work,
    }
    for channel, error in measured.items():
        assert error < limits[channel], f"{channel}: {error:.6g} >= {limits[channel]}"

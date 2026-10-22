"""Lock a direct MATLAB heading-bank branch split without promoting the case."""

import gzip
import hashlib
import json
from pathlib import Path

import numpy as np

from wecsim import PMWave, WEC, WorldPoint


FIXTURES = Path(__file__).parent / "fixtures/matlab"
HYDRO = FIXTURES / "oswec_multi_wave/oswec.h5"
BANK = np.arange(-30, 30.25, .25)


def test_direct_matlab_80s_initial_yaw_branch_split():
    """Only the 1e-7 rad initial yaw changes; both runs select their own bank."""
    with np.load(FIXTURES / "oswec_yaw_80_initial_perturb_source.npz") as saved:
        source = saved["trajectory"]
        source_wave = saved["wave"]
    with np.load(FIXTURES / "oswec_yaw_160_same_bank_half.npz") as saved:
        baseline = saved["half"][:len(source)]
        baseline_wave = saved["half_wave"][:len(source)]
    evidence = json.loads((FIXTURES / "oswec_yaw_80_initial_perturb_evidence.json").read_text())
    with gzip.open(FIXTURES / "oswec_yaw_80_initial_perturb_input.m.gz", "rb") as file:
        assert hashlib.sha256(file.read()).hexdigest() == evidence["input_sha256"]
    assert source.shape == (16_001, 7)
    assert source_wave.shape == baseline_wave.shape == (16_001, 2)
    assert evidence["published_case_promotion"] is False
    np.testing.assert_allclose(source_wave, baseline_wave, rtol=0, atol=1e-12)
    source_mismatch = np.flatnonzero(source[:, 4] != baseline[:, 4])
    assert len(source_mismatch) == evidence["heading_mismatch_samples"] == 7
    assert source[source_mismatch[0], 0] == evidence["first_heading_mismatch_s"] == 68.33

    wec = WEC("OSWEC selected-bank initial-yaw source diagnostic")
    flap = wec.body("flap", HYDRO, mass=12_700, inertia=(1.85e6,) * 3,
                    passive_yaw=True, yaw_heading_bank=BANK)
    wec.body("base", HYDRO, mass=999, inertia=(999,) * 3)
    yaw = wec.coordinate("yaw", flap.move("yaw", pivot=WorldPoint(0, 0, -8.9)))
    wec.rotational_pto("hinge", yaw, damping=120_000)
    result = wec.run(
        PMWave(2.5, 8, direction=10, seed=1, phase_generator="matlab"),
        dt=.005, end_time=80, ramp_time=100, radiation_memory=40,
        solver="wecsim_ode4", initial_coordinate={"yaw": 1e-7},
    )
    np.testing.assert_allclose(result.time, source[:, 0], rtol=0, atol=1e-10)
    np.testing.assert_allclose(result.wave_elevation, source_wave[:, 1],
                               rtol=0, atol=1e-12)
    angle = result.bodies["flap"].position[:, 5]
    speed = result.bodies["flap"].velocity[:, 5]
    relative = 10 - np.degrees(angle)
    chosen = BANK[np.argmin(np.abs(BANK[:, None] - relative[None, :]), axis=0)]
    np.testing.assert_array_equal(chosen, source[:, 4])
    assert np.max(np.abs(angle - source[:, 1])) < 1e-12
    assert np.max(np.abs(speed - source[:, 2])) < 1e-12
    assert np.max(np.abs(result.ptos["hinge"].force - source[:, 5])) < 1e-6
    source_work = -np.trapezoid(source[:, 6], source[:, 0])
    python_work = np.trapezoid(result.ptos["hinge"].absorbed_power, result.time)
    assert source_work > 0
    assert abs(python_work - source_work) / source_work < 1e-10

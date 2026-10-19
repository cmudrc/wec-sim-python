"""Keep the full-duration OWC source-step separation visible in CI.

The Python trajectory advances independently. This diagnostic does not turn
the source's step-sensitive orifice case into a physical-motion parity gate.
"""

from pathlib import Path

import numpy as np

from wecsim import OrificePTO, PMWave, WEC


FIXTURES = Path(__file__).parent / "fixtures/matlab"
SOURCE = FIXTURES / "owc_orifice_full_maxstep/source.npz"
QUARTER = FIXTURES / "owc_orifice_full_maxstep/quarter_source.npz"
EIGHTH = FIXTURES / "owc_orifice_full_maxstep/eighth_source.npz"
PUBLISHED = FIXTURES / "owc_orifice_130s/source.npz"
HYDRO = FIXTURES / "owc_orifice.h5"


def test_full_owc_source_step_audit(tmp_path):
    with np.load(SOURCE, allow_pickle=False) as saved, \
            np.load(QUARTER, allow_pickle=False) as quarter_saved, \
            np.load(EIGHTH, allow_pickle=False) as eighth_saved, \
            np.load(PUBLISHED, allow_pickle=False) as published:
        trace = {key: saved[key] for key in saved.files}
        quarter = {key: quarter_saved[key] for key in quarter_saved.files}
        eighth = {key: eighth_saved[key] for key in eighth_saved.files}
        phases = published["phases"].copy()
        parameters = published["orifice_parameters"].copy()
    time = trace["time"]
    assert time.shape == (26_001,)
    assert trace["nominal_components"].shape == (500, 4)
    np.testing.assert_allclose(time, np.arange(26_001) * .005,
                               rtol=0, atol=1e-10)
    np.testing.assert_array_equal(trace["nominal_components"],
                                  trace["half_components"])
    np.testing.assert_array_equal(trace["nominal_wave"], trace["half_wave"])
    np.testing.assert_allclose(trace["nominal_components"][:, 3], phases,
                               rtol=0, atol=1e-14)
    assert all(trace[f"{arm}_motion"].shape == (26_001, 4)
               and trace[f"{arm}_power"].shape == time.shape
               for arm in ("nominal", "half"))
    np.testing.assert_array_equal(quarter["time"], time)
    np.testing.assert_array_equal(quarter["components"],
                                  trace["nominal_components"])
    np.testing.assert_array_equal(quarter["wave"], trace["nominal_wave"])
    assert float(quarter["simu_dt"]) == .005
    assert float(quarter["max_step"]) == .00125
    assert quarter["motion"].shape == quarter["velocity"].shape == (26_001, 4)
    assert np.isfinite(quarter["motion"]).all()
    assert int(quarter["orifice_samples"]) > 400_000
    np.testing.assert_array_equal(eighth["time"], time)
    np.testing.assert_array_equal(eighth["components"],
                                  quarter["components"])
    np.testing.assert_array_equal(eighth["wave"], quarter["wave"])
    assert eighth["motion"].shape == (26_001, 4)
    assert np.isfinite(eighth["motion"]).all()
    assert float(eighth["simu_dt"]) == .005
    assert float(eighth["max_step"]) == .000625
    assert int(eighth["orifice_samples"]) > int(quarter["orifice_samples"])
    with np.load(FIXTURES / "owc_orifice_dt_30s/maxstep_only.npz",
                 allow_pickle=False) as prefix:
        np.testing.assert_array_equal(quarter["motion"][:6001],
                                      prefix["motion"])
        np.testing.assert_array_equal(quarter["velocity"][:6001],
                                      prefix["velocity"])

    phase_file = tmp_path / "phases.csv"
    np.savetxt(phase_file, phases, delimiter=",")
    wec = WEC("OWC full-duration source-step diagnostic")
    body = wec.body("OWC", HYDRO, inertia=(99.28, 11.04, 99.2))
    wec.floating_gbm(
        body, orifice=OrificePTO(*parameters),
        orifice_force_path="published_owc",
        heave_linear_damping=100, mode_linear_damping=100,
        heave_drag_cd=1.2, heave_drag_area=8,
        pitch_drag_cd=1.2, pitch_drag_area=8,
    )
    wave = PMWave(height=1, period=4, phase_file=phase_file)
    result = wec.run(
        wave,
        dt=.005, end_time=130, ramp_time=10, radiation_memory=15,
    )
    np.testing.assert_allclose(result.time, time, rtol=0, atol=1e-10)
    np.testing.assert_allclose(result.wave_elevation, trace["nominal_wave"],
                               rtol=0, atol=1e-11)

    python_motion = np.column_stack((
        result.bodies["OWC"].position[:, 0],
        result.bodies["OWC"].position[:, 2],
        result.coordinates["pitch_unwrapped"].position,
        result.flexible_modes["OWC"].position[:, 0],
    ))
    nominal = trace["nominal_motion"]
    half = trace["half_motion"]
    self_change = np.max(np.abs(nominal - half), axis=0)
    nominal_gap = np.sqrt(np.mean((python_motion - nominal)**2, axis=0))
    half_gap = np.sqrt(np.mean((python_motion - half)**2, axis=0))
    print("OWC MATLAB MaxStep self-change:", self_change)
    print("OWC Python/source RMS gaps, nominal and half:",
          nominal_gap, half_gap)

    # The immutable source traces show step sensitivity over all 130 s.
    # Bound Python's current residual without requiring it to remain large:
    # a future improvement should pass this test, too.
    assert self_change[0] > .03
    assert self_change[2] > .004
    assert half_gap[2] < .002
    assert half_gap[0] < .14

    # One more source MaxStep halving changes full-duration surge substantially;
    # the smaller Python gap is a diagnostic, not a published-motion gate.
    quarter_motion = quarter["motion"]
    quarter_source_change = np.max(np.abs(quarter_motion - half), axis=0)
    quarter_gap = np.abs(python_motion - quarter_motion)
    quarter_rms = np.sqrt(np.mean(quarter_gap**2, axis=0))
    print("OWC MATLAB quarter MaxStep self-change:", quarter_source_change)
    print("OWC Python/quarter source peak and RMS gaps:",
          np.max(quarter_gap, axis=0), quarter_rms)
    assert quarter_source_change[0] > .12
    assert np.max(quarter_gap[:, 0]) < .08
    assert quarter_rms[0] < .06

    # The next source step moves surge another 66 mm. Its smaller RMS gap
    # to Python is still a diagnostic, not a converged published-case gate.
    eighth_motion = eighth["motion"]
    eighth_source_change = np.max(np.abs(eighth_motion - quarter_motion), axis=0)
    eighth_gap = np.abs(python_motion - eighth_motion)
    eighth_rms = np.sqrt(np.mean(eighth_gap**2, axis=0))
    print("OWC MATLAB eighth MaxStep self-change:", eighth_source_change)
    print("OWC Python/eighth source peak and RMS gaps:",
          np.max(eighth_gap, axis=0), eighth_rms)
    assert .06 < eighth_source_change[0] < .08
    assert np.max(eighth_gap[:, 0]) < .11
    assert eighth_rms[0] < .04

    python_work = np.trapezoid(result.ptos["orifice"].absorbed_power, time)
    nominal_work = np.trapezoid(trace["nominal_power"], time)
    half_work = np.trapezoid(trace["half_power"], time)
    assert min(python_work, nominal_work, half_work) > 0
    source_work_change = abs(half_work - nominal_work) / nominal_work
    half_work_gap = abs(python_work - half_work) / half_work
    print("OWC source work self-change and Python/half gap:",
          source_work_change, half_work_gap)
    assert source_work_change < .01
    assert half_work_gap < .06
    quarter_work = float(quarter["work_internal"])
    quarter_source_work_change = abs(quarter_work - half_work) / half_work
    quarter_work_gap = abs(python_work - quarter_work) / quarter_work
    print("OWC quarter source work self-change and Python gap:",
          quarter_source_work_change, quarter_work_gap)
    assert quarter_source_work_change < .005
    assert quarter_work_gap < .05
    eighth_work = float(eighth["work_internal"])
    eighth_source_work_change = abs(eighth_work - quarter_work) / quarter_work
    eighth_work_gap = abs(python_work - eighth_work) / eighth_work
    print("OWC eighth source work self-change and Python gap:",
          eighth_source_work_change, eighth_work_gap)
    assert .002 < eighth_source_work_change < .003
    assert eighth_work_gap < .05

    # Refine only the Python solver. Its shared-sample trajectory and work
    # should be stable even though the MATLAB/Python case gap remains open.
    fine = wec.run(
        wave, dt=.0025, end_time=130, ramp_time=10, radiation_memory=15,
    )
    np.testing.assert_allclose(fine.time[::2], time, rtol=0, atol=1e-10)
    fine_motion = np.column_stack((
        fine.bodies["OWC"].position[::2, 0],
        fine.bodies["OWC"].position[::2, 2],
        fine.coordinates["pitch_unwrapped"].position[::2],
        fine.flexible_modes["OWC"].position[::2, 0],
    ))
    python_step_change = np.max(np.abs(python_motion - fine_motion), axis=0)
    print("OWC Python half-step self-change:", python_step_change)
    assert np.all(python_step_change < [.0001, .002, .00002, .0003])
    fine_work = np.trapezoid(fine.ptos["orifice"].absorbed_power, fine.time)
    assert abs(fine_work - python_work) / python_work < .001

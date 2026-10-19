"""Keep the full-duration OWC source-step separation visible in CI.

The Python trajectory advances independently. This diagnostic does not turn
the source's step-sensitive orifice case into a physical-motion parity gate.
"""

from pathlib import Path

import numpy as np

from wecsim import OrificePTO, PMWave, WEC


FIXTURES = Path(__file__).parent / "fixtures/matlab"
SOURCE = FIXTURES / "owc_orifice_full_maxstep/source.npz"
PUBLISHED = FIXTURES / "owc_orifice_130s/source.npz"
HYDRO = FIXTURES / "owc_orifice.h5"


def test_full_owc_step_change_does_not_close_source_route_gap(tmp_path):
    with np.load(SOURCE, allow_pickle=False) as saved, \
            np.load(PUBLISHED, allow_pickle=False) as published:
        trace = {key: saved[key] for key in saved.files}
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
    result = wec.run(
        PMWave(height=1, period=4, phase_file=phase_file),
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

    # The source's 10 s step sensitivity persists over 130 s. Refining its
    # solver materially narrows pitch, but not the remaining surge gap.
    assert self_change[0] > .03
    assert self_change[2] > .004
    assert half_gap[2] < .6 * nominal_gap[2]
    assert half_gap[0] > .1

    python_work = np.trapezoid(result.ptos["orifice"].absorbed_power, time)
    nominal_work = np.trapezoid(trace["nominal_power"], time)
    half_work = np.trapezoid(trace["half_power"], time)
    assert min(python_work, nominal_work, half_work) > 0
    source_work_change = abs(half_work - nominal_work) / nominal_work
    half_work_gap = abs(python_work - half_work) / half_work
    print("OWC source work self-change and Python/half gap:",
          source_work_change, half_work_gap)
    assert source_work_change < .01
    assert half_work_gap > .04

"""Pair the physical OWC path with a one-connection corrected MATLAB model.

The published Simulink application omits piston force from its flexible
state integrator. This derived comparison does not count as published-case
parity; it checks the physical Python path against a corrected source model.
"""

from pathlib import Path

import numpy as np

from wecsim import OrificePTO, PMWave, WEC


FIXTURES = Path(__file__).parent / "fixtures/matlab"
SOURCE = FIXTURES / "owc_orifice_corrected_piston/source.npz"
PUBLISHED = FIXTURES / "owc_orifice_130s/source.npz"
ORIGINAL_FINE = FIXTURES / "owc_orifice_full_maxstep/sixteenth_source.npz"
HYDRO = FIXTURES / "owc_orifice.h5"
MOTION_LIMITS = np.array([.020, .010, .001, .002])


def test_full_duration_corrected_source_and_physical_python(tmp_path):
    with (np.load(SOURCE, allow_pickle=False) as corrected,
          np.load(PUBLISHED, allow_pickle=False) as published,
          np.load(ORIGINAL_FINE, allow_pickle=False) as original):
        source = {key: corrected[key] for key in corrected.files}
        phases = published["phases"].copy()
        parameters = published["orifice_parameters"].copy()
        original_motion = original["motion"].copy()
        np.testing.assert_array_equal(source["time"], original["time"])
        np.testing.assert_array_equal(source["wave"], original["wave"])
        np.testing.assert_array_equal(source["components"],
                                      original["components"])
        assert source["components"].shape == (500, 4)
        np.testing.assert_allclose(source["components"][:, 3], phases,
                                   rtol=0, atol=1e-14)

    time = source["time"]
    np.testing.assert_allclose(time, np.arange(26_001) * .005,
                               rtol=0, atol=1e-10)
    assert float(source["coarse_max_step"]) == .0003125
    assert float(source["fine_max_step"]) == .00015625
    assert int(source["coarse_orifice_samples"]) == 416_001
    assert int(source["fine_orifice_samples"]) == 832_001
    assert float(source["coarse_max_mach"]) < float(parameters[4])
    assert float(source["fine_max_mach"]) < float(parameters[4])
    assert float(source["coarse_flag_fraction"]) == 0
    assert float(source["fine_flag_fraction"]) == 0
    assert source["coarse_motion"].shape == (26_001, 4)
    assert source["fine_motion"].shape == (26_001, 4)
    assert np.isfinite(source["coarse_motion"]).all()
    assert np.isfinite(source["fine_motion"]).all()

    # These full-duration screens were set before the refined corrected-source
    # runs. A close Python trace alone is insufficient if MATLAB is unstable.
    source_change = np.max(np.abs(source["fine_motion"]
                                  - source["coarse_motion"]), axis=0)
    coarse_work = float(source["coarse_work_internal"])
    fine_work = float(source["fine_work_internal"])
    assert coarse_work > 0 and fine_work > 0
    source_work_change = abs(fine_work - coarse_work) / coarse_work
    print("corrected OWC MATLAB self-change:", source_change,
          "PTO work:", source_work_change)
    assert np.all(source_change <= MOTION_LIMITS), source_change
    assert source_work_change <= .005

    phase_file = tmp_path / "phases.csv"
    np.savetxt(phase_file, phases, delimiter=",")
    wec = WEC("OWC physical piston-path comparison")
    body = wec.body("OWC", HYDRO, inertia=(99.28, 11.04, 99.2))
    wec.floating_gbm(
        body, orifice=OrificePTO(*parameters),
        heave_linear_damping=100, mode_linear_damping=100,
        heave_drag_cd=1.2, heave_drag_area=8,
        pitch_drag_cd=1.2, pitch_drag_area=8,
    )
    result = wec.run(
        PMWave(height=1, period=4, phase_file=phase_file),
        dt=.005, end_time=130, ramp_time=10, radiation_memory=15,
    )
    np.testing.assert_allclose(result.time, time, rtol=0, atol=1e-10)
    np.testing.assert_allclose(result.wave_elevation, source["wave"],
                               rtol=0, atol=1e-11)
    python_motion = np.column_stack((
        result.bodies["OWC"].position[:, 0],
        result.bodies["OWC"].position[:, 2],
        result.coordinates["pitch_unwrapped"].position,
        result.flexible_modes["OWC"].position[:, 0],
    ))
    python_gap = np.max(np.abs(python_motion - source["fine_motion"]),
                        axis=0)
    python_work = np.trapezoid(result.ptos["orifice"].absorbed_power, time)
    python_work_gap = abs(python_work - fine_work) / fine_work
    print("physical Python/corrected MATLAB peak gap:", python_gap,
          "PTO work:", python_work_gap)
    assert np.all(python_gap <= MOTION_LIMITS), python_gap
    assert python_work_gap <= .01

    # The original and corrected models share the sea and solver step. Their
    # large separation isolates the omitted piston contribution in the source.
    original_change = np.max(np.abs(
        source["coarse_motion"] - original_motion), axis=0)
    assert original_change[1] > 1
    assert original_change[3] > .1

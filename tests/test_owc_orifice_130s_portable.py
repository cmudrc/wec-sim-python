"""Keep the published 130 s OWC source-convention comparison in ordinary CI."""

from pathlib import Path

import numpy as np

from wecsim import OrificePTO, PMWave, WEC


FIXTURE = Path(__file__).parent / "fixtures/matlab/owc_orifice_130s/source.npz"
HYDRO = Path(__file__).parent / "fixtures/matlab/owc_orifice.h5"


def _within(actual, source, limit, label):
    assert actual.shape == source.shape, label
    assert np.isfinite(actual).all() and np.isfinite(source).all(), label
    gap = float(np.max(np.abs(actual - source)))
    assert gap < limit, f"{label}: {gap:.6g} exceeds {limit}"


def test_published_owc_source_path_full_duration(tmp_path):
    # The archived traces are outputs only. Python builds its own wave and
    # advances rigid motion, flexible motion, and the PTO from the HDF5 data.
    with np.load(FIXTURE, allow_pickle=False) as saved:
        source = {key: saved[key] for key in saved.files}
    time = source["time"]
    assert time.shape == (26_001,)
    np.testing.assert_allclose(time, np.arange(len(time)) * 0.005,
                               rtol=0, atol=1e-10)
    assert source["phases"].shape == (500,)

    phase_file = tmp_path / "phases.csv"
    np.savetxt(phase_file, source["phases"], delimiter=",")
    wec = WEC("OWC published force path")
    body = wec.body("OWC", HYDRO, inertia=(99.28, 11.04, 99.2))
    wec.floating_gbm(
        body, orifice=OrificePTO(*source["orifice_parameters"]),
        orifice_force_path="published_owc",
        heave_linear_damping=100, mode_linear_damping=100,
        heave_drag_cd=1.2, heave_drag_area=8,
        pitch_drag_cd=1.2, pitch_drag_area=8,
    )
    result = wec.run(
        PMWave(height=1, period=4, phase_file=phase_file),
        dt=0.005, end_time=130, ramp_time=10, radiation_memory=15,
    )
    assert result.case["constraint"]["orifice_force_path"] == "published_owc"
    np.testing.assert_allclose(result.time, time, rtol=0, atol=1e-10)
    _within(result.wave_elevation, source["wave_elevation"], 1e-11,
            "incident sea")

    pose = result.bodies["OWC"].position
    velocity = result.bodies["OWC"].velocity
    flex = result.flexible_modes["OWC"]
    _within(pose[:, 0], source["body_position"][:, 0], 0.175, "surge")
    _within(pose[:, 2], source["body_position"][:, 1], 0.038, "heave")
    _within(result.coordinates["pitch_unwrapped"].position,
            source["body_position"][:, 2], 0.012, "unwrapped pitch")
    _within(velocity[:, 0], source["body_velocity"][:, 0], 0.008,
            "surge speed")
    _within(velocity[:, 2], source["body_velocity"][:, 1], 0.080,
            "heave speed")
    _within(velocity[:, 4], source["body_velocity"][:, 2], 0.0045,
            "pitch speed")
    _within(flex.position[:, 0], source["flex_position"], 0.006,
            "flexible displacement")
    _within(flex.velocity[:, 0], source["flex_velocity"], 0.051,
            "flexible speed")

    pto = result.ptos["orifice"]
    source_force = source["orifice_force"]
    source_power = source["orifice_power"]
    _within(pto.force, source_force,
            0.06 * np.max(np.abs(source_force)), "PTO reaction")
    _within(pto.absorbed_power, source_power,
            0.08 * np.max(np.abs(source_power)), "PTO absorbed power")
    source_work = np.trapezoid(source_power, time)
    model_work = np.trapezoid(pto.absorbed_power, time)
    assert source_work > 0
    assert abs(model_work - source_work) < 0.06 * source_work

    first_ten = time <= 10
    _within(pose[first_ten, 2], source["body_position"][first_ten, 1],
            0.002, "first-ten heave")
    _within(flex.position[first_ten, 0], source["flex_position"][first_ten],
            0.0005, "first-ten flexible displacement")

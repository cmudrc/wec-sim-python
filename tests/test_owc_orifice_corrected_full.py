"""Gate the independently advanced OWC against the corrected source model.

The corrected model adds piston force to the flexible state input. Its source
trace is derived from the published application, so this is not a published
input badge gate.
"""

import os
from pathlib import Path

import numpy as np

from wecsim import OrificePTO, PMWave, WEC


FIXTURES = Path(__file__).parent / "fixtures/matlab"
SOURCE = Path(os.environ.get(
    "WEC_SIM_OWC_CORRECTED_SOURCE",
    FIXTURES / "owc_orifice_corrected_piston/full_source.npz",
))
ORIGINAL = FIXTURES / "owc_orifice_corrected_piston/source.npz"
PUBLISHED = FIXTURES / "owc_orifice_130s/source.npz"
HYDRO = FIXTURES / "owc_orifice.h5"


def _peak(actual, expected, limit, name):
    assert actual.shape == expected.shape, name
    assert np.isfinite(actual).all() and np.isfinite(expected).all(), name
    error = float(np.max(np.abs(actual - expected)))
    assert error <= limit, f"{name}: {error:.8g} > {limit:.8g}"


def test_corrected_owc_full_motion_and_air_system(tmp_path):
    with np.load(SOURCE, allow_pickle=False) as saved:
        source = {name: saved[name] for name in saved.files}
    with np.load(ORIGINAL, allow_pickle=False) as saved:
        existing = {name: saved[name] for name in saved.files}
    with np.load(PUBLISHED, allow_pickle=False) as saved:
        phases = saved["phases"].copy()
        parameters = saved["orifice_parameters"].copy()

    body = source["body"]
    flex = source["flex"]
    air = source["orifice"]
    time = body[:, 0]
    _peak(time, existing["time"], 1e-10, "output time")
    _peak(source["wave"][:, 1], existing["wave"], 1e-10,
          "incident wave")
    _peak(source["components"], existing["components"], 1e-12,
          "PM components and phases")
    _peak(source["components"][:, 3], phases, 1e-12, "PM phase input")
    assert not np.any(air[:, 2])
    assert source["summary"][1] < parameters[4]

    # The exported fluid signals must agree with the source piston speed and
    # orifice law before comparing independently advanced Python outputs.
    orifice = OrificePTO(*parameters)
    source_air = orifice.evaluate(flex[:, 2])
    _peak(source_air.force, air[:, 1], 1e-5, "source PTO force law")
    _peak(source_air.pressure_drop / 1000, air[:, 3], 1e-7,
          "source pressure law")
    _peak(source_air.flow_rate, air[:, 4], 1e-10, "source flow law")
    _peak(source_air.absorbed_power / 1000, air[:, 5], 1e-7,
          "source power law")

    phase_file = tmp_path / "phases.csv"
    np.savetxt(phase_file, phases, delimiter=",")
    wec = WEC("OWC corrected piston-path gate")
    floating = wec.body("OWC", HYDRO, inertia=(99.28, 11.04, 99.2))
    wec.floating_gbm(
        floating, orifice=orifice,
        heave_linear_damping=100, mode_linear_damping=100,
        heave_drag_cd=1.2, heave_drag_area=8,
        pitch_drag_cd=1.2, pitch_drag_area=8,
    )
    result = wec.run(
        PMWave(height=1, period=4, phase_file=phase_file),
        dt=.005, end_time=130, ramp_time=10, radiation_memory=15,
    )
    _peak(result.time, time, 1e-10, "Python time")
    _peak(result.wave_elevation, source["wave"][:, 1], 1e-10,
          "Python incident wave")

    pose = result.bodies["OWC"].position
    velocity = result.bodies["OWC"].velocity
    flexible = result.flexible_modes["OWC"]
    source_pitch = np.r_[0, np.cumsum(
        (body[1:, 11] + body[:-1, 11]) * .005 / 2,
    )]
    source_motion = np.column_stack((
        body[:, 1], body[:, 3], source_pitch, flex[:, 1],
    ))
    source_reproduction = np.max(np.abs(
        source_motion - existing["fine_motion"]), axis=0)
    assert np.all(source_reproduction <= [.020, .010, .001, .002]), (
        source_reproduction
    )
    assert abs(float(source["summary"][0])
               - float(existing["fine_work_internal"])) / float(
                   existing["fine_work_internal"]) <= .005
    _peak(pose[:, 0], body[:, 1], .020, "surge")
    _peak(pose[:, 2], body[:, 3], .010, "heave")
    _peak(result.coordinates["pitch_unwrapped"].position,
          source_pitch, .001, "pitch")
    _peak(flexible.position[:, 0], flex[:, 1], .002,
          "flexible displacement")
    _peak(velocity[:, 0], body[:, 7], .020, "surge speed")
    _peak(velocity[:, 2], body[:, 9], .020, "heave speed")
    _peak(velocity[:, 4], body[:, 11], .002, "pitch speed")
    _peak(flexible.velocity[:, 0], flex[:, 2], .005, "piston speed")

    pto = result.ptos["orifice"]
    python_air = orifice.evaluate(pto.velocity)
    force_scale = np.max(np.abs(air[:, 1]))
    pressure_scale = np.max(np.abs(air[:, 3]))
    flow_scale = np.max(np.abs(air[:, 4]))
    power_scale = np.max(np.abs(air[:, 5]))
    assert min(force_scale, pressure_scale, flow_scale, power_scale) > 0
    _peak(pto.force, air[:, 1], .03 * force_scale, "PTO force")
    _peak(python_air.pressure_drop / 1000, air[:, 3],
          .03 * pressure_scale, "air pressure")
    _peak(python_air.flow_rate, air[:, 4], .03 * flow_scale, "air flow")
    _peak(pto.absorbed_power / 1000, air[:, 5],
          .04 * power_scale, "PTO power")
    assert not np.any(python_air.compressibility_flag)
    source_work = float(source["summary"][0])
    python_work = float(np.trapezoid(pto.absorbed_power, time))
    assert abs(python_work - source_work) / source_work <= .01

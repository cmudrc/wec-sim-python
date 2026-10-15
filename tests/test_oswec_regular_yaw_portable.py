"""Gate the published 600 s OSWEC selected-heading regular-wave case."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim import RegularWave, WEC, WorldPoint
from wecsim.hydrodynamics.bodyClass import BodyClass
from wecsim.hydrodynamics.passiveYaw import (
    NearestHeadingExcitation, PassiveYawExcitation,
)


ROOT = Path(__file__).resolve().parents[1]
PORTABLE = Path(__file__).parent / "fixtures/matlab/oswec_regular_yaw/source.npz"
APPLICATIONS = os.environ.get("WEC_SIM_APPLICATIONS_DIR")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")


def _source():
    if APPLICATIONS and REFERENCE:
        source = Path(REFERENCE)
        prefix = "OSWEC_VARIABLE_YAW_PUBLISHED"
        flap = np.loadtxt(source / f"{prefix}_regular_body1.csv", delimiter=",")
        base = np.loadtxt(source / f"{prefix}_regular_body2.csv", delimiter=",")
        pto = np.loadtxt(source / f"{prefix}_regular_pto1.csv", delimiter=",")
        wave = np.loadtxt(source / f"{prefix}_wave.csv", delimiter=",")
        selected = np.loadtxt(source / f"{prefix}_selected_heading.csv", delimiter=",")
        assert flap.shape == base.shape == pto.shape == (60001, 25)
        assert wave.shape == (60001, 2)
        assert selected.shape == (60001, 3)
        for other in (base[:, 0], pto[:, 0], wave[:, 0], selected[:, 0]):
            np.testing.assert_allclose(other, flap[:, 0], rtol=0, atol=1e-10)
        values = {
            "time": flap[:, 0], "yaw": flap[:, 6],
            "yaw_speed": flap[:, 12], "excitation": flap[:, 19:25],
            "base_position": base[:, 1:7], "base_velocity": base[:, 7:13],
            "pto_force": pto[:, 17], "pto_velocity": pto[:, 11],
            "pto_power": pto[:, 23], "wave_elevation": wave[:, 1],
            "selected_heading": selected[:, 2],
        }
        hydro = (Path(APPLICATIONS)
                 / "_Common_Input_Files/OSWEC/hydroData/oswec.h5")
        return values, hydro
    if not APPLICATIONS and not REFERENCE:
        with np.load(PORTABLE) as saved:
            values = {name: saved[name] for name in saved.files}
        return values, ROOT / "tests/fixtures/matlab/oswec_multi_wave/oswec.h5"
    pytest.fail("provide both WEC_SIM_APPLICATIONS_DIR and "
                "WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")


def test_published_oswec_regular_selected_yaw_against_matlab():
    source, hydro = _source()
    assert source["time"].shape == (60001,)
    assert source["excitation"].shape == (60001, 6)
    assert all(np.isfinite(value).all() for value in source.values())
    headings = np.arange(-30, 30.25, .25)

    wec = WEC("OSWEC regular variable yaw")
    flap = wec.body(
        "flap", hydro, mass=12700, inertia=(1.85e6,) * 3,
        passive_yaw=True, yaw_heading_bank=headings,
    )
    wec.body("base", hydro, mass=999, inertia=(999,) * 3)
    yaw = wec.coordinate("yaw", flap.move("yaw", pivot=WorldPoint(0, 0, -8.9)))
    wec.rotational_pto("hinge", yaw, damping=120000)
    result = wec.run(RegularWave(2.5, 8, 10), dt=.01,
                     end_time=600, ramp_time=100)

    np.testing.assert_allclose(result.time, source["time"], rtol=0, atol=1e-10)
    np.testing.assert_allclose(result.wave_elevation, source["wave_elevation"],
                               rtol=0, atol=1e-12)
    np.testing.assert_allclose(result.bodies["base"].position,
                               source["base_position"], rtol=0, atol=1e-10)
    np.testing.assert_allclose(result.bodies["base"].velocity,
                               source["base_velocity"], rtol=0, atol=1e-10)

    body = BodyClass(str(hydro))
    body.bodyNumber = 1
    body.bodyTotal = 2
    body.readH5file()
    model = PassiveYawExcitation.from_hydro_data(
        body.hydroData, omega=2 * np.pi / 8,
        incident_direction=10, amplitude=1.25, ramp_time=100,
        rho=1000, g=9.81, spline_frequency=True,
    )
    bank = NearestHeadingExcitation(model, headings)
    source_headings = np.array([bank.heading(angle) for angle in source["yaw"]])
    np.testing.assert_array_equal(source_headings, source["selected_heading"])
    python_headings = np.array([
        bank.heading(angle) for angle in result.bodies["flap"].position[:, 5]
    ])
    assert np.count_nonzero(np.diff(source["selected_heading"])) == 374
    np.testing.assert_array_equal(python_headings, source["selected_heading"])
    source_path_force = np.stack([
        bank.force(t, angle) for t, angle in zip(result.time, source["yaw"])
    ])
    assert np.max(np.abs(source_path_force - source["excitation"])) < 1e-6

    motion = result.bodies["flap"]
    assert np.max(np.abs(motion.position[:, 5] - source["yaw"])) < 1e-6
    assert np.max(np.abs(motion.velocity[:, 5] - source["yaw_speed"])) < 1e-6
    assert np.max(np.abs(result.ptos["hinge"].force - source["pto_force"])) < .01
    np.testing.assert_allclose(source["pto_force"],
                               -120000 * source["pto_velocity"],
                               rtol=0, atol=1e-5)
    source_work = -np.trapezoid(source["pto_power"], source["time"])
    python_work = np.trapezoid(result.ptos["hinge"].absorbed_power,
                               result.time)
    assert source_work > 0
    assert abs(python_work - source_work) / source_work < 1e-6

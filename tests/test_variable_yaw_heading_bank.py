"""Nearest-heading selection used by the published variable-hydro yaw case."""

import os
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pytest

import wecsim.caseDynamics as case_dynamics
from wecsim import PMWave, RegularWave, WEC, WorldPoint
from wecsim.irregularWave import pm_equal_energy_components
from wecsim.passiveYaw import (
    NearestHeadingExcitation, NearestSampledHeadingExcitation,
    PassiveYawExcitation, SampledPassiveYawExcitation,
)


def _model():
    headings = np.arange(0.0, 360.0, 10.0)
    real = np.zeros((6, len(headings)))
    real[0] = 10
    real[1] = 20
    real[5] = headings
    return PassiveYawExcitation(
        headings, real, np.zeros_like(real),
        incident_direction=10, omega=1, amplitude=1, ramp_time=0,
    )


def test_nearest_bank_uses_source_relative_angle_and_first_tie():
    bank = NearestHeadingExcitation(_model(), np.arange(-40, 41, 2))
    assert bank.heading(0) == 10
    assert bank.heading(np.deg2rad(9)) == 0
    assert bank.heading(np.deg2rad(-31)) == 40
    assert bank.heading(np.deg2rad(51)) == -40
    assert bank.heading(np.deg2rad(9)) == 0  # a one-degree tie chooses lower

    yaw = np.deg2rad(9.2)
    actual = bank.force(0, yaw)
    expected = _model().force(0, 0, coefficient_heading=0)
    np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-13)
    assert actual[5] == 0


def test_nearest_bank_rejects_invalid_grid():
    model = _model()
    for headings in ([0], [0, 0], [2, 1], [-181, 0], [0, np.nan]):
        with pytest.raises(ValueError, match="heading bank"):
            NearestHeadingExcitation(model, headings)


def test_nearest_sampled_bank_selects_heading_without_yaw_rotation():
    directions = np.arange(0, 360, 10)
    force_grid = np.zeros((2, 1, len(directions), 6))
    force_grid[:, 0, 0, 0] = 10
    force_grid[:, 0, 0, 1] = 20
    force_grid[:, 0, 0, 5] = 30
    sampled = SampledPassiveYawExcitation(
        np.array([0, .01]), directions, np.array([10]),
        force_grid, np.zeros(2),
    )
    bank = NearestSampledHeadingExcitation(sampled, np.arange(-40, 41, 2))
    yaw = np.deg2rad(9.2)
    assert bank.heading(yaw) == 0
    np.testing.assert_allclose(
        bank.force(0, yaw), [10, 20, 0, 0, 0, 30], rtol=0, atol=1e-13,
    )


@pytest.mark.skipif(
    not (os.environ.get("WEC_SIM_APPLICATIONS_DIR")
         and os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")),
    reason="fresh variable-yaw MATLAB output not provided",
)
def test_two_degree_bank_against_pinned_matlab():
    _check_bank_against_pinned_matlab(
        "OSWEC_VARIABLE_YAW_2DEG", "regular_2deg", np.arange(-40, 41, 2),
    )


@pytest.mark.skipif(
    not (os.environ.get("WEC_SIM_APPLICATIONS_DIR")
         and os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")),
    reason="fresh variable-yaw MATLAB output not provided",
)
def test_published_bank_against_pinned_matlab():
    _check_bank_against_pinned_matlab(
        "OSWEC_VARIABLE_YAW_PUBLISHED", "regular",
        np.arange(-30, 30.25, .25),
    )


@pytest.mark.skipif(
    not (os.environ.get("WEC_SIM_APPLICATIONS_DIR")
         and os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")),
    reason="fresh variable-yaw MATLAB output not provided",
)
def test_irregular_120s_bank_against_pinned_matlab():
    """Pair the source force path and independent 120 s PM trajectory."""
    applications = Path(os.environ["WEC_SIM_APPLICATIONS_DIR"])
    source = Path(os.environ["WEC_SIM_MATLAB_MODEL_OUTPUT_DIR"])
    hydro = applications / "_Common_Input_Files/OSWEC/hydroData/oswec.h5"
    prefix = "OSWEC_VARIABLE_YAW_IRREGULAR_120S"
    flap = np.loadtxt(source / f"{prefix}_irregular_120s_body1.csv", delimiter=",")
    base = np.loadtxt(source / f"{prefix}_irregular_120s_body2.csv", delimiter=",")
    pto = np.loadtxt(source / f"{prefix}_irregular_120s_pto1.csv", delimiter=",")
    wave = np.loadtxt(source / f"{prefix}_wave.csv", delimiter=",")
    selected = np.loadtxt(source / f"{prefix}_selected_heading.csv", delimiter=",")
    realization = np.loadtxt(source / f"{prefix}_components.csv", delimiter=",")
    assert flap.shape == base.shape == pto.shape == (12001, 25)
    assert wave.shape == (12001, 2)
    assert selected.shape == (12001, 3)
    assert realization.shape == (500, 4)

    components = pm_equal_energy_components(
        hydro, significant_height=2.5, peak_period=8,
        directions=[10], spreading=[1], seed=1, phase_generator="matlab",
    )
    for actual, expected in (
        (components.omega, realization[:, 0]),
        (components.spectral_amplitude, realization[:, 1]),
        (components.d_omega, realization[:, 2]),
        (components.phase[:, 0], realization[:, 3]),
    ):
        np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-12)

    sampled = SampledPassiveYawExcitation.from_hydro_data(
        _hydro_data(hydro), components, dt=.01, end_time=120,
        ramp_time=100, rho=1000, g=9.81,
    )
    headings = np.arange(-30, 30.25, .25)
    bank = NearestSampledHeadingExcitation(sampled, headings)
    np.testing.assert_allclose(sampled.elevation, wave[:, 1], rtol=0, atol=1e-12)
    source_heading = np.array([bank.heading(angle) for angle in flap[:, 6]])
    np.testing.assert_array_equal(source_heading, selected[:, 2])
    source_path_force = np.stack([
        bank.force(t, angle) for t, angle in zip(flap[:, 0], flap[:, 6])
    ])
    assert np.max(np.abs(source_path_force - flap[:, 19:25])) < 1e-6

    wec = WEC("OSWEC irregular variable yaw")
    moving = wec.body(
        "flap", hydro, mass=12700, inertia=(1.85e6,) * 3,
        passive_yaw=True, yaw_heading_bank=headings,
    )
    wec.body("base", hydro, mass=999, inertia=(999,) * 3)
    yaw = wec.coordinate("yaw", moving.move("yaw", pivot=WorldPoint(0, 0, -8.9)))
    wec.rotational_pto("hinge", yaw, damping=120000)
    result = wec.run(
        PMWave(2.5, 8, direction=10, seed=1, phase_generator="matlab"),
        dt=.01, end_time=120, ramp_time=100, radiation_memory=40,
    )
    np.testing.assert_allclose(result.time, flap[:, 0], rtol=0, atol=1e-10)
    np.testing.assert_allclose(result.wave_elevation, wave[:, 1],
                               rtol=0, atol=1e-12)
    np.testing.assert_allclose(result.bodies["base"].position,
                               base[:, 1:7], rtol=0, atol=1e-10)
    np.testing.assert_allclose(result.bodies["base"].velocity,
                               base[:, 7:13], rtol=0, atol=1e-10)
    assert np.max(np.abs(result.bodies["flap"].position[:, 5] - flap[:, 6])) < .002
    assert np.max(np.abs(result.bodies["flap"].velocity[:, 5] - flap[:, 12])) < .0003
    assert np.max(np.abs(result.ptos["hinge"].force - pto[:, 17])) < 35
    np.testing.assert_allclose(pto[:, 17], -120000 * pto[:, 11],
                               rtol=0, atol=1e-5)
    source_work = -np.trapezoid(pto[:, 23], flap[:, 0])
    python_work = np.trapezoid(result.ptos["hinge"].absorbed_power, result.time)
    assert source_work > 0
    assert abs(python_work - source_work) / source_work < .005


@pytest.mark.skipif(
    not os.environ.get("WEC_SIM_APPLICATIONS_DIR"),
    reason="pinned variable-yaw HDF5 not provided",
)
def test_continuous_irregular_120s_step_stability():
    """Check a continuous-heading physical control on the pinned PM sea."""
    hydro = (Path(os.environ["WEC_SIM_APPLICATIONS_DIR"])
             / "_Common_Input_Files/OSWEC/hydroData/oswec.h5")
    results = []
    for dt in (.01, .005):
        wec = WEC("OSWEC continuous-heading stability control")
        flap = wec.body(
            "flap", hydro, mass=12700, inertia=(1.85e6,) * 3,
            passive_yaw=True,
        )
        wec.body("base", hydro, mass=999, inertia=(999,) * 3)
        yaw = wec.coordinate(
            "yaw", flap.move("yaw", pivot=WorldPoint(0, 0, -8.9)),
        )
        wec.rotational_pto("hinge", yaw, damping=120000)
        results.append(wec.run(
            PMWave(2.5, 8, direction=10, seed=1,
                   phase_generator="matlab"),
            dt=dt, end_time=120, ramp_time=100, radiation_memory=40,
        ))

    coarse, fine = results
    assert coarse.time.shape == (12001,)
    assert fine.time.shape == (24001,)
    np.testing.assert_allclose(coarse.time, fine.time[::2], rtol=0, atol=1e-10)
    np.testing.assert_allclose(coarse.wave_elevation,
                               fine.wave_elevation[::2], rtol=0, atol=1e-12)
    for kind, limit in (("position", 2e-5), ("velocity", 3e-6)):
        a = getattr(coarse.bodies["flap"], kind)[:, 5]
        b = getattr(fine.bodies["flap"], kind)[::2, 5]
        assert np.max(np.abs(a - b)) < limit, kind
    work = [np.trapezoid(result.ptos["hinge"].absorbed_power, result.time)
            for result in results]
    assert min(work) > 0
    assert abs(work[0] - work[1]) / work[1] < 1e-4


@pytest.mark.skipif(
    not (os.environ.get("WEC_SIM_APPLICATIONS_DIR")
         and os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")),
    reason="fresh variable-yaw MATLAB output not provided",
)
def test_irregular_600s_source_force_against_pinned_matlab():
    """Check the published full-length sea and force path on source motion."""
    applications = Path(os.environ["WEC_SIM_APPLICATIONS_DIR"])
    source = Path(os.environ["WEC_SIM_MATLAB_MODEL_OUTPUT_DIR"])
    hydro = applications / "_Common_Input_Files/OSWEC/hydroData/oswec.h5"
    prefix = "OSWEC_VARIABLE_YAW_IRREGULAR_600S"
    flap = np.loadtxt(source / f"{prefix}_irregular_600s_body1.csv", delimiter=",")
    pto = np.loadtxt(source / f"{prefix}_irregular_600s_pto1.csv", delimiter=",")
    wave = np.loadtxt(source / f"{prefix}_wave.csv", delimiter=",")
    selected = np.loadtxt(source / f"{prefix}_selected_heading.csv", delimiter=",")
    realization = np.loadtxt(source / f"{prefix}_components.csv", delimiter=",")
    assert flap.shape == pto.shape == (60001, 25)
    assert wave.shape == (60001, 2)
    assert selected.shape == (60001, 3)
    assert realization.shape == (500, 4)

    components = pm_equal_energy_components(
        hydro, significant_height=2.5, peak_period=8,
        directions=[10], spreading=[1], seed=1, phase_generator="matlab",
    )
    for actual, expected in (
        (components.omega, realization[:, 0]),
        (components.spectral_amplitude, realization[:, 1]),
        (components.d_omega, realization[:, 2]),
        (components.phase[:, 0], realization[:, 3]),
    ):
        np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-12)

    sampled = SampledPassiveYawExcitation.from_hydro_data(
        _hydro_data(hydro), components, dt=.01, end_time=600,
        ramp_time=100, rho=1000, g=9.81,
    )
    bank = NearestSampledHeadingExcitation(sampled, np.arange(-30, 30.25, .25))
    np.testing.assert_allclose(sampled.time, flap[:, 0], rtol=0, atol=1e-10)
    np.testing.assert_allclose(sampled.elevation, wave[:, 1], rtol=0, atol=1e-12)
    source_heading = np.array([bank.heading(angle) for angle in flap[:, 6]])
    np.testing.assert_array_equal(source_heading, selected[:, 2])
    source_path_force = np.stack([
        bank.force(t, angle) for t, angle in zip(flap[:, 0], flap[:, 6])
    ])
    assert np.max(np.abs(source_path_force - flap[:, 19:25])) < 1e-6
    np.testing.assert_allclose(pto[:, 17], -120000 * pto[:, 11],
                               rtol=0, atol=1e-5)


@pytest.mark.skipif(
    not (os.environ.get("WEC_SIM_APPLICATIONS_DIR")
         and os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")),
    reason="fresh variable-yaw MATLAB output not provided",
)
def test_irregular_600s_prescribed_heading_dynamics_diagnostic():
    """Advance Python motion with only the source's saved heading events."""
    applications = Path(os.environ["WEC_SIM_APPLICATIONS_DIR"])
    source = Path(os.environ["WEC_SIM_MATLAB_MODEL_OUTPUT_DIR"])
    hydro = applications / "_Common_Input_Files/OSWEC/hydroData/oswec.h5"
    prefix = "OSWEC_VARIABLE_YAW_IRREGULAR_600S"
    flap = np.loadtxt(source / f"{prefix}_irregular_600s_body1.csv", delimiter=",")
    pto = np.loadtxt(source / f"{prefix}_irregular_600s_pto1.csv", delimiter=",")
    selected = np.loadtxt(source / f"{prefix}_selected_heading.csv", delimiter=",")
    headings = np.arange(-30, 30.25, .25)
    assert flap.shape == pto.shape == (60001, 25)
    assert selected.shape == (60001, 3)

    class SourceHeadingSchedule:
        def __init__(self, model, bank_headings):
            self.model = model
            np.testing.assert_array_equal(bank_headings, headings)

        def force(self, at_time, yaw):
            index = min(max(round(at_time / .01), 0), 60000)
            coefficient_yaw = np.deg2rad(
                self.model.incident_directions[0] - selected[index, 2]
            )
            return self.model.force(at_time, 0.0,
                                    coefficient_yaw=coefficient_yaw)

    wec = WEC("OSWEC variable yaw: prescribed heading diagnostic")
    moving = wec.body(
        "flap", hydro, mass=12700, inertia=(1.85e6,) * 3,
        passive_yaw=True, yaw_heading_bank=headings,
    )
    wec.body("base", hydro, mass=999, inertia=(999,) * 3)
    yaw = wec.coordinate("yaw", moving.move("yaw", pivot=WorldPoint(0, 0, -8.9)))
    wec.rotational_pto("hinge", yaw, damping=120000)
    with patch.object(case_dynamics, "NearestSampledHeadingExcitation",
                      SourceHeadingSchedule):
        result = wec.run(
            PMWave(2.5, 8, direction=10, seed=1,
                   phase_generator="matlab"),
            dt=.01, end_time=600, ramp_time=100, radiation_memory=40,
        )

    np.testing.assert_allclose(result.time, flap[:, 0], rtol=0, atol=1e-10)
    assert np.max(np.abs(result.bodies["flap"].position[:, 5] - flap[:, 6])) < .01
    assert np.max(np.abs(result.bodies["flap"].velocity[:, 5] - flap[:, 12])) < .0001
    assert np.max(np.abs(result.ptos["hinge"].force - pto[:, 17])) < 10
    source_work = -np.trapezoid(pto[:, 23], flap[:, 0])
    replay_work = np.trapezoid(result.ptos["hinge"].absorbed_power, result.time)
    assert source_work > 0
    assert abs(replay_work - source_work) / source_work < .0002


def _check_bank_against_pinned_matlab(prefix, case, headings):
    applications = Path(os.environ["WEC_SIM_APPLICATIONS_DIR"])
    source = Path(os.environ["WEC_SIM_MATLAB_MODEL_OUTPUT_DIR"])
    flap = np.loadtxt(source / f"{prefix}_{case}_body1.csv", delimiter=",")
    base = np.loadtxt(source / f"{prefix}_{case}_body2.csv", delimiter=",")
    pto = np.loadtxt(source / f"{prefix}_{case}_pto1.csv", delimiter=",")
    selected = np.loadtxt(source / f"{prefix}_selected_heading.csv", delimiter=",")
    wave = np.loadtxt(source / f"{prefix}_wave.csv", delimiter=",")
    hydro = applications / "_Common_Input_Files/OSWEC/hydroData/oswec.h5"

    wec = WEC("OSWEC variable yaw")
    moving = wec.body(
        "flap", hydro, mass=12700, inertia=(1.85e6,) * 3,
        passive_yaw=True, yaw_heading_bank=headings,
    )
    wec.body("base", hydro, mass=999, inertia=(999,) * 3)
    yaw = wec.coordinate("yaw", moving.move("yaw", pivot=WorldPoint(0, 0, -8.9)))
    wec.rotational_pto("hinge", yaw, damping=120000)
    result = wec.run(RegularWave(2.5, 8, 10), dt=0.01,
                     end_time=600, ramp_time=100)

    assert flap.shape == base.shape == pto.shape == (60001, 25)
    assert selected.shape == (60001, 3)
    assert wave.shape == (60001, 2)
    np.testing.assert_allclose(result.time, flap[:, 0], rtol=0, atol=1e-10)
    np.testing.assert_allclose(result.wave_elevation, wave[:, 1], rtol=0, atol=1e-12)
    np.testing.assert_allclose(result.bodies["base"].position,
                               base[:, 1:7], rtol=0, atol=1e-10)
    np.testing.assert_allclose(result.bodies["base"].velocity,
                               base[:, 7:13], rtol=0, atol=1e-10)

    model = PassiveYawExcitation.from_hydro_data(
        _hydro_data(hydro), omega=2 * np.pi / 8,
        incident_direction=10, amplitude=1.25, ramp_time=100,
        rho=1000, g=9.81, spline_frequency=True,
    )
    bank = NearestHeadingExcitation(model, headings)
    source_heading = np.array([bank.heading(angle) for angle in flap[:, 6]])
    np.testing.assert_array_equal(source_heading, selected[:, 2])

    source_path_force = np.stack([
        bank.force(t, angle) for t, angle in zip(result.time, flap[:, 6])
    ])
    assert np.max(np.abs(source_path_force - flap[:, 19:25])) < 1e-6
    assert np.max(np.abs(result.bodies["flap"].position[:, 5]
                         - flap[:, 6])) < 1e-6
    assert np.max(np.abs(result.bodies["flap"].velocity[:, 5]
                         - flap[:, 12])) < 1e-6
    assert np.max(np.abs(result.ptos["hinge"].force
                         - pto[:, 17])) < .01
    source_work = -np.trapezoid(pto[:, 23], flap[:, 0])
    python_work = np.trapezoid(result.ptos["hinge"].absorbed_power, result.time)
    assert abs(python_work - source_work) / abs(source_work) < 1e-6


def _hydro_data(path):
    from wecsim.bodyClass import BodyClass

    body = BodyClass(str(path))
    body.bodyNumber = 1
    body.bodyTotal = 2
    body.readH5file()
    return body.hydroData

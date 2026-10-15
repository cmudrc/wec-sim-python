"""Published held-yaw source numerics through the public WEC API in ordinary CI."""

from pathlib import Path

import numpy as np
import pytest

from wecsim import PMWave, WEC, WorldPoint
from wecsim.environment.irregularWave import pm_equal_energy_components
from wecsim.hydrodynamics.bodyClass import BodyClass
from wecsim.hydrodynamics.passiveYaw import (
    HeldPassiveYawExcitation, SampledPassiveYawExcitation,
)


FIXTURES = Path(__file__).resolve().parent / "fixtures/matlab"
HYDRO = FIXTURES / "oswec_multi_wave/oswec.h5"
SOURCE = FIXTURES / "oswec_held_yaw/source.npz"
HEADING_EVENTS = {1: 141, 2: 52, 3: 91}


def test_default_held_yaw_commits_heading_and_reports_consistent_work():
    """The physical implicit path must commit each accepted heading once."""
    wec = WEC("OSWEC implicit held yaw")
    flap = wec.body(
        "flap", HYDRO, mass=12700, inertia=(1.85e6,) * 3,
        passive_yaw=True, passive_yaw_threshold=1,
    )
    wec.body("base", HYDRO, mass=999, inertia=(999,) * 3)
    yaw = wec.coordinate("yaw", flap.move("yaw", pivot=WorldPoint(0, 0, -8.9)))
    wec.rotational_pto("hinge", yaw, damping=120000)
    result = wec.run(
        PMWave(2.5, 8, direction=10, seed=1, phase_generator="matlab"),
        dt=.05, end_time=2, ramp_time=0, radiation_memory=1,
        initial_coordinate={"yaw": .02},
    )
    assert result.time.shape == (41,)
    assert np.ptp(result.coordinates["yaw"].position) > .02

    # Replay only the accepted Python yaw states through the held-heading
    # force model. This checks force logging and commit order, not MATLAB
    # trajectory parity for the step-sensitive published case.
    body = BodyClass(str(HYDRO))
    body.bodyNumber = 1
    body.bodyTotal = 2
    body.readH5file()
    components = pm_equal_energy_components(
        HYDRO, significant_height=2.5, peak_period=8,
        directions=[10], spreading=[1], count=500,
        seed=1, phase_generator="matlab",
    )
    sampled = SampledPassiveYawExcitation.from_hydro_data(
        body.hydroData, components, dt=.05, end_time=2,
        ramp_time=0, rho=1000, g=9.81,
    )
    held = HeldPassiveYawExcitation(sampled, threshold=1)
    headings = []
    for instant, angle in zip(result.time,
                              result.coordinates["yaw"].position, strict=True):
        held.commit(instant, np.array([angle]))
        headings.append(held.last_heading)
    assert np.count_nonzero(np.abs(np.diff(headings)) > 1e-6) == 1
    np.testing.assert_allclose(
        dict(result.raw.extra_outputs)["body1_excitation_force"],
        held.force_history, rtol=0, atol=1e-8,
    )
    pto = result.ptos["hinge"]
    np.testing.assert_allclose(pto.force, -120000 * pto.velocity,
                               rtol=0, atol=1e-10)
    np.testing.assert_allclose(pto.absorbed_power,
                               120000 * pto.velocity**2,
                               rtol=0, atol=1e-10)
    assert np.trapezoid(pto.absorbed_power, result.time) > 0


def test_implicit_and_source_ode4_approach_the_same_short_yaw_motion():
    """Refining the step should remove their numerical, not physical, gap.

    The long published held-yaw run is sensitive to heading-event timing.
    This short seeded window checks the two independent solvers at the same
    resolution before that accumulated event error dominates.
    """
    gaps = []
    for dt in (.05, .0125, .003125):
        results = {}
        for solver in ("implicit", "wecsim_ode4"):
            wec = WEC("OSWEC short held-yaw solver convergence")
            flap = wec.body(
                "flap", HYDRO, mass=12700, inertia=(1.85e6,) * 3,
                passive_yaw=True, passive_yaw_threshold=1,
            )
            wec.body("base", HYDRO, mass=999, inertia=(999,) * 3)
            yaw = wec.coordinate(
                "yaw", flap.move("yaw", pivot=WorldPoint(0, 0, -8.9)),
            )
            wec.rotational_pto("hinge", yaw, damping=120000)
            results[solver] = wec.run(
                PMWave(2.5, 8, direction=10, seed=1,
                       phase_generator="matlab"),
                dt=dt, end_time=2, ramp_time=0, radiation_memory=1,
                initial_coordinate={"yaw": .02}, solver=solver,
            )
        implicit, source = results["implicit"], results["wecsim_ode4"]
        np.testing.assert_array_equal(implicit.time, source.time)
        gaps.append((
            np.max(np.abs(implicit.coordinates["yaw"].position
                          - source.coordinates["yaw"].position)),
            np.max(np.abs(implicit.coordinates["yaw"].velocity
                          - source.coordinates["yaw"].velocity)),
        ))
    gaps = np.asarray(gaps)
    assert np.isfinite(gaps).all()
    # Each fourfold step refinement must reduce both motion gaps by at least
    # a factor of three; the finest run also has an absolute motion bound.
    assert np.all(gaps[1:] < gaps[:-1] / 3), gaps
    assert gaps[-1, 0] < 1e-6
    assert gaps[-1, 1] < 1.5e-6


def _heading_events(hydro_data, components, time, yaw):
    sampled = SampledPassiveYawExcitation.from_hydro_data(
        hydro_data, components, dt=.01, end_time=250,
        ramp_time=100, rho=1000, g=9.81,
    )
    held = HeldPassiveYawExcitation(sampled, threshold=1)
    headings = np.empty(len(time))
    for index, instant in enumerate(time):
        held.commit(instant, np.array([yaw[index]]))
        headings[index] = held.last_heading
    return np.flatnonzero(np.abs(np.diff(headings)) > 1e-5)


@pytest.mark.parametrize("seed", (1, 2, 3))
def test_published_250s_held_yaw_source_path_from_phase_seed(seed):
    """The opt-in solver reproduces source output; the physical default is separate."""
    with np.load(SOURCE) as archive:
        time = archive["time"]
        source = {key.removeprefix(f"seed{seed}_"): archive[key]
                  for key in archive.files if key.startswith(f"seed{seed}_")}
    components = pm_equal_energy_components(
        HYDRO, significant_height=2.5, peak_period=8,
        directions=[10], spreading=[1], count=500,
        seed=seed, phase_generator="matlab",
    )
    generated = np.column_stack((
        components.omega, components.spectral_amplitude,
        components.d_omega, components.phase[:, 0],
    ))
    np.testing.assert_allclose(generated, source["components"],
                               rtol=0, atol=1e-12)

    wec = WEC("OSWEC published held yaw")
    flap = wec.body("flap", HYDRO, mass=12700, inertia=(1.85e6,) * 3,
                    passive_yaw=True, passive_yaw_threshold=1)
    wec.body("base", HYDRO, mass=999, inertia=(999,) * 3)
    yaw = wec.coordinate("yaw", flap.move("yaw", pivot=WorldPoint(0, 0, -8.9)))
    wec.rotational_pto("hinge", yaw, damping=120000)
    result = wec.run(
        PMWave(2.5, 8, direction=10, seed=seed, phase_generator="matlab"),
        dt=.01, end_time=250, ramp_time=100, radiation_memory=40,
        solver="wecsim_ode4",
    )

    np.testing.assert_allclose(result.time, time, rtol=0, atol=1e-10)
    np.testing.assert_allclose(result.wave_elevation,
                               source["wave_elevation"], rtol=0, atol=1e-11)
    expected_base = np.array([0, 0, -10.9, 0, 0, 0])
    np.testing.assert_allclose(result.bodies["base"].position,
                               np.broadcast_to(expected_base, (len(time), 6)),
                               rtol=0, atol=1e-12)
    np.testing.assert_allclose(result.bodies["base"].velocity,
                               np.zeros((len(time), 6)), rtol=0, atol=1e-12)
    actual_force = dict(result.raw.extra_outputs)["body1_excitation_force"]
    assert np.max(np.abs(actual_force[:, 5] - source["yaw_excitation"])) < .001

    angle = result.bodies["flap"].position[:, 5]
    speed = result.bodies["flap"].velocity[:, 5]
    pto = result.ptos["hinge"]
    assert np.max(np.abs(angle - source["yaw"])) < 1e-6
    assert np.max(np.abs(speed - source["yaw_speed"])) < 1e-7
    assert np.max(np.abs(pto.force - source["pto_torque"])) < .01
    np.testing.assert_allclose(source["pto_torque"],
                               -120000 * source["pto_speed"],
                               rtol=0, atol=1e-5)
    source_work = -np.trapezoid(source["pto_power"], time)
    python_work = np.trapezoid(pto.absorbed_power, time)
    assert source_work > 0
    assert abs(python_work - source_work) / source_work < 1e-6

    body = BodyClass(str(HYDRO))
    body.bodyNumber = 1
    body.bodyTotal = 2
    body.readH5file()
    source_events = _heading_events(body.hydroData, components,
                                    time, source["yaw"])
    python_events = _heading_events(body.hydroData, components, time, angle)
    assert len(source_events) == HEADING_EVENTS[seed]
    np.testing.assert_array_equal(python_events, source_events)

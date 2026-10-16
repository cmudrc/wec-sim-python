"""Measure ode45 MaxStep independently of End_Stops sampling settings."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim import LinearHardStops, RegularWave, WEC


OUTPUT = os.environ.get("WEC_SIM_RM3_END_STOPS_MAXSTEP_DIR")
pytestmark = pytest.mark.skipif(not OUTPUT, reason="MATLAB audit output absent")
END_TIME = int(os.environ.get("WEC_SIM_RM3_END_STOPS_MAXSTEP_END_TIME", "120"))
PINNED = (Path(__file__).parent / "fixtures/matlab"
          / "rm3_end_stops_published_prefix.npz")
SMALLER_120 = (Path(__file__).parent / "fixtures/matlab"
               / "rm3_end_stops_maxstep_120s.npz")
REFINED = (Path(__file__).parent / "fixtures/matlab"
           / "rm3_end_stops_full.npz")


def _read(directory, label):
    count = round(END_TIME / .1) + 1
    settings = np.loadtxt(directory / f"{label}_settings.csv", delimiter=",")
    wave = np.loadtxt(directory / f"{label}_wave.csv", delimiter=",")
    pto = np.loadtxt(directory / f"{label}_pto.csv", delimiter=",")
    bodies = [np.loadtxt(directory / f"{label}_body{body}.csv",
                         delimiter=",") for body in (1, 2)]
    assert settings.shape == (4,)
    assert wave.shape == (count, 2)
    assert pto.shape == (count, 4)
    assert all(body.shape == (count, 7) for body in bodies)
    assert all(np.isfinite(array).all() for array in (wave, pto, *bodies))
    return settings, wave, pto, bodies


def test_maxstep_only_source_separation():
    """Pin the unchanged source, then report whether MaxStep closes its gap."""
    directory = Path(OUTPUT)
    published = _read(directory, "published")
    smaller = _read(directory, "maxstep_0025")
    assert END_TIME in (120, 400)
    np.testing.assert_allclose(published[0], (.1, END_TIME, 100, .1),
                               rtol=0, atol=1e-12)
    np.testing.assert_allclose(smaller[0], (.1, END_TIME, 100, .025),
                               rtol=0, atol=1e-12)
    time = np.arange(round(END_TIME / .1) + 1) * .1
    for wave, pto, bodies in ((published[1], published[2], published[3]),
                              (smaller[1], smaller[2], smaller[3])):
        np.testing.assert_allclose(wave[:, 0], time, rtol=0, atol=1e-9)
        np.testing.assert_allclose(pto[:, 0], time, rtol=0, atol=1e-9)
        for body in bodies:
            np.testing.assert_allclose(body[:, 0], time, rtol=0, atol=1e-9)
    np.testing.assert_allclose(published[1], smaller[1],
                               rtol=0, atol=1e-12)

    with np.load(PINNED, allow_pickle=False) as pinned:
        np.testing.assert_allclose(pinned["time"], time[:1201],
                                   rtol=0, atol=1e-10)
        for column, limit in enumerate((1e-7, 1e-6, 1e-2)):
            np.testing.assert_allclose(
                published[2][:1201, column + 1], pinned["pto"][:, column],
                rtol=0, atol=limit,
            )
        for index, body in enumerate(published[3]):
            np.testing.assert_allclose(
                body[:1201, 1:4], pinned["body_position"][:, index],
                rtol=0, atol=1e-7,
            )
            np.testing.assert_allclose(
                body[:1201, 4:7], pinned["body_velocity"][:, index],
                rtol=0, atol=1e-6,
            )
    if END_TIME == 400:
        # The changed full-duration run must retain the already gated 120 s
        # prefix before a late contact event is interpreted as a new effect.
        with np.load(SMALLER_120, allow_pickle=False) as prefix:
            for column, limit in enumerate((1e-7, 1e-6, 1e-2)):
                np.testing.assert_allclose(
                    smaller[2][:1201, column + 1],
                    prefix["pto"][:, column], rtol=0, atol=limit,
                )
            for index, body in enumerate(smaller[3]):
                np.testing.assert_allclose(
                    body[:1201, 1:4], prefix["body_position"][:, index],
                    rtol=0, atol=1e-7,
                )
                np.testing.assert_allclose(
                    body[:1201, 4:7], prefix["body_velocity"][:, index],
                    rtol=0, atol=1e-6,
                )

    with np.load(REFINED, allow_pickle=False) as refined:
        fine_time = refined["dt00125_time"][::8][:len(time)]
        fine_pto = refined["dt00125_pto"][::8][:len(time)]
    np.testing.assert_allclose(fine_time, time, rtol=0, atol=1e-10)
    after_contact = time >= 68.0
    errors = {}
    for label, arm in (("published", published),
                       ("maxstep_0025", smaller)):
        pto = arm[2]
        error = pto[:, 1] - fine_pto[:, 0]
        rms = np.sqrt(np.mean(error[after_contact] ** 2))
        maximum = np.max(np.abs(error[after_contact]))
        damper_work = np.trapezoid(1_200_000 * pto[:, 2] ** 2, time)
        fine_work = np.trapezoid(1_200_000 * fine_pto[:, 1] ** 2, time)
        work_gap = abs(damper_work - fine_work) / fine_work
        errors[label] = (rms, maximum, work_gap)
        print(f"{label}: post-contact stroke RMS={rms:.9g} m, "
              f"maximum={maximum:.9g} m, damper work gap="
              f"{(damper_work - fine_work) / fine_work:.6%}")

    # Predeclared source self-change gate: the one-setting edit must remove
    # most of the published-to-refined discrepancy through the whole window.
    for improved, original, fraction in zip(
            errors["maxstep_0025"], errors["published"], (.1, .1, .25),
            strict=True):
        assert original > 0
        assert improved < fraction * original

    # This audit tests a solver mechanism. Its result does not qualify the
    # published case or alter the physical Python default on its own.
    assert np.count_nonzero(np.abs(published[2][:, 1]) > .6) > 0
    assert np.count_nonzero(np.abs(smaller[2][:, 1]) > .6) > 0


@pytest.mark.skipif(
    END_TIME != 400 or not os.environ.get("WEC_SIM_RM3_H5"),
    reason="full-duration MATLAB output and pinned RM3 HDF5 required",
)
def test_independent_public_motion_against_full_maxstep_source():
    """Advance the published 400 s inputs without source states or loads."""
    source = _read(Path(OUTPUT), "maxstep_0025")
    time = source[2][:, 0]
    hydro = Path(os.environ["WEC_SIM_RM3_H5"])
    wec = WEC("RM3 End Stops full-duration MaxStep audit")
    float_body = wec.body("float", hydro,
                          inertia=(0, 21_306_090.66, 0))
    spar = wec.body("spar", hydro,
                    inertia=(0, 94_407_091.24, 0))
    wec.floating_joint(
        float_body, spar, damping=1_200_000,
        hard_stops=LinearHardStops(-.6, .6, 1e8, 1e8),
    )
    result = wec.run(RegularWave(2.5, 8), dt=.1, end_time=400,
                     ramp_time=100)
    np.testing.assert_allclose(result.time, time, rtol=0, atol=1e-9)
    np.testing.assert_allclose(result.wave_elevation, source[1][:, 1],
                               rtol=0, atol=1e-8)

    for index, name in enumerate(("float", "spar")):
        body = result.bodies[name]
        for actual, expected, limits, field in (
            (body.position, source[3][index][:, 1:4],
             (.005, .006, 1e-4), "position"),
            (body.velocity, source[3][index][:, 4:7],
             (.005, .03, 1e-4), "velocity"),
        ):
            for column, (axis, limit) in enumerate(zip(
                    (0, 2, 4), limits, strict=True)):
                error = np.max(np.abs(actual[:, axis] - expected[:, column]))
                print(f"{name} {field} axis {axis}: {error:.9g}")
                assert error < limit

    pto = result.ptos["relative_heave"]
    for actual, expected, limit, label in (
        (pto.stroke, source[2][:, 1], .006, "PTO stroke"),
        (pto.velocity, source[2][:, 2], .035, "PTO speed"),
        (pto.force, source[2][:, 3], 450_000, "PTO force"),
    ):
        error = np.max(np.abs(actual - expected))
        print(f"{label}: {error:.9g}")
        assert error < limit
    python_work = np.trapezoid(pto.absorbed_power, result.time)
    source_work = np.trapezoid(
        1_200_000 * source[2][:, 2] ** 2, time,
    )
    assert source_work > 0
    work_gap = abs(python_work - source_work) / source_work
    print(f"relative damper work: {work_gap:.9g}")
    assert work_gap < .007

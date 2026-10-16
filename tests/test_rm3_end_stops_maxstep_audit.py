"""Measure ode45 MaxStep independently of End_Stops sampling settings."""

import os
from pathlib import Path

import numpy as np
import pytest


OUTPUT = os.environ.get("WEC_SIM_RM3_END_STOPS_MAXSTEP_DIR")
pytestmark = pytest.mark.skipif(not OUTPUT, reason="MATLAB audit output absent")
PINNED = (Path(__file__).parent / "fixtures/matlab"
          / "rm3_end_stops_published_prefix.npz")
REFINED = (Path(__file__).parent / "fixtures/matlab"
           / "rm3_end_stops_full.npz")


def _read(directory, label):
    settings = np.loadtxt(directory / f"{label}_settings.csv", delimiter=",")
    wave = np.loadtxt(directory / f"{label}_wave.csv", delimiter=",")
    pto = np.loadtxt(directory / f"{label}_pto.csv", delimiter=",")
    bodies = [np.loadtxt(directory / f"{label}_body{body}.csv",
                         delimiter=",") for body in (1, 2)]
    assert settings.shape == (4,)
    assert wave.shape == (1201, 2)
    assert pto.shape == (1201, 4)
    assert all(body.shape == (1201, 7) for body in bodies)
    assert all(np.isfinite(array).all() for array in (wave, pto, *bodies))
    return settings, wave, pto, bodies


def test_maxstep_only_source_separation():
    """Pin the unchanged source, then report whether MaxStep closes its gap."""
    directory = Path(OUTPUT)
    published = _read(directory, "published")
    smaller = _read(directory, "maxstep_0025")
    np.testing.assert_allclose(published[0], (.1, 120, 100, .1),
                               rtol=0, atol=1e-12)
    np.testing.assert_allclose(smaller[0], (.1, 120, 100, .025),
                               rtol=0, atol=1e-12)
    time = np.arange(1201) * .1
    for wave, pto, bodies in ((published[1], published[2], published[3]),
                              (smaller[1], smaller[2], smaller[3])):
        np.testing.assert_allclose(wave[:, 0], time, rtol=0, atol=1e-9)
        np.testing.assert_allclose(pto[:, 0], time, rtol=0, atol=1e-9)
        for body in bodies:
            np.testing.assert_allclose(body[:, 0], time, rtol=0, atol=1e-9)
    np.testing.assert_allclose(published[1], smaller[1],
                               rtol=0, atol=1e-12)

    with np.load(PINNED, allow_pickle=False) as pinned:
        np.testing.assert_allclose(pinned["time"], time, rtol=0, atol=1e-10)
        for column, limit in enumerate((1e-7, 1e-6, 1e-2)):
            np.testing.assert_allclose(
                published[2][:, column + 1], pinned["pto"][:, column],
                rtol=0, atol=limit,
            )
        for index, body in enumerate(published[3]):
            np.testing.assert_allclose(
                body[:, 1:4], pinned["body_position"][:, index],
                rtol=0, atol=1e-7,
            )
            np.testing.assert_allclose(
                body[:, 4:7], pinned["body_velocity"][:, index],
                rtol=0, atol=1e-6,
            )

    with np.load(REFINED, allow_pickle=False) as refined:
        fine_time = refined["dt00125_time"][::8][:1201]
        fine_pto = refined["dt00125_pto"][::8][:1201]
    np.testing.assert_allclose(fine_time, time, rtol=0, atol=1e-10)
    after_contact = time >= 68.0
    for label, arm in (("published", published),
                       ("maxstep_0025", smaller)):
        pto = arm[2]
        error = pto[:, 1] - fine_pto[:, 0]
        rms = np.sqrt(np.mean(error[after_contact] ** 2))
        maximum = np.max(np.abs(error[after_contact]))
        damper_work = np.trapezoid(1_200_000 * pto[:, 2] ** 2, time)
        fine_work = np.trapezoid(1_200_000 * fine_pto[:, 1] ** 2, time)
        print(f"{label}: post-contact stroke RMS={rms:.9g} m, "
              f"maximum={maximum:.9g} m, damper work gap="
              f"{(damper_work - fine_work) / fine_work:.6%}")

    # This audit tests a solver mechanism. Its result does not qualify the
    # published case or alter the physical Python default on its own.
    assert np.count_nonzero(np.abs(published[2][:, 1]) > .6) > 0
    assert np.count_nonzero(np.abs(smaller[2][:, 1]) > .6) > 0

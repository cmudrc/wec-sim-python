"""An RM3 device's motion must not depend on its coordinate basis."""

from pathlib import Path

import numpy as np
from scipy.integrate import trapezoid

from wecsim import NoWave, RegularWave, WEC


RM3 = Path(__file__).resolve().parents[1] / "examples/data/rm3.h5"


def _run_rm3(*, mixed_basis: bool, body_to_body: bool = True,
             regular_wave: bool = False):
    wec = WEC("RM3 heave basis", body_to_body=body_to_body)
    float_body = wec.body("float", RM3)
    spar = wec.body("spar", RM3)
    if mixed_basis:
        # Physical heaves are common + relative and common - relative.
        # Both coordinates therefore move both bodies, including through
        # the cross-body added mass and radiation blocks in the RM3 HDF5.
        wec.coordinate("common", float_body.move("heave"), spar.move("heave"))
        wec.coordinate("relative", float_body.move("heave"),
                       spar.move("heave", scale=-1))
        position = {"common": .1, "relative": .3}
        speed = {"common": .025, "relative": .075}
    else:
        wec.coordinate("float_heave", float_body.move("heave"))
        wec.coordinate("spar_heave", spar.move("heave"))
        position = {"float_heave": .4, "spar_heave": -.2}
        speed = {"float_heave": .1, "spar_heave": -.05}
    wec.pto("main", float_body.at(0, 0, 0), spar.at(0, 0, 0),
            axis=(0, 0, 1), stiffness=25_000, damping=100_000)
    wave = RegularWave(2.5, 8) if regular_wave else NoWave()
    simulation = ({"dt": .05, "end_time": 2., "ramp_time": 0.}
                  if regular_wave else
                  {"dt": .01, "end_time": .5, "radiation_memory": .2})
    return wec.run(wave, **simulation,
                   initial_coordinate=position, initial_speed=speed)


def test_two_body_radiation_and_pto_work_are_coordinate_basis_invariant():
    direct = _run_rm3(mixed_basis=False)
    mixed = _run_rm3(mixed_basis=True)

    for name in ("float", "spar"):
        for field in ("position", "velocity"):
            np.testing.assert_allclose(
                getattr(mixed.bodies[name], field),
                getattr(direct.bodies[name], field), rtol=0, atol=1e-11,
            )
    for field in ("stroke", "velocity", "force", "absorbed_power"):
        np.testing.assert_allclose(
            getattr(mixed.ptos["main"], field),
            getattr(direct.ptos["main"], field), rtol=0, atol=1e-8,
        )
    assert np.ptp(direct.ptos["main"].stroke) > .01
    work = trapezoid(direct.ptos["main"].absorbed_power, x=direct.time)
    assert work > 0
    np.testing.assert_allclose(
        trapezoid(mixed.ptos["main"].absorbed_power, x=mixed.time),
        work, rtol=0, atol=1e-8,
    )

    # The identity check must exercise actual cross-body hydrodynamics.
    uncoupled = _run_rm3(mixed_basis=False, body_to_body=False)
    assert np.max(np.abs(
        direct.bodies["spar"].position
        - uncoupled.bodies["spar"].position,
    )) > 1e-4


def test_regular_wave_cross_body_damping_is_coordinate_basis_invariant():
    """The same physical RM3 motion must survive a change of heave basis.

    This exercises the regular-wave cross-body radiation damping matrix,
    rather than the convolution kernel used by the no-wave test above.
    """
    direct = _run_rm3(mixed_basis=False, regular_wave=True)
    mixed = _run_rm3(mixed_basis=True, regular_wave=True)
    uncoupled = _run_rm3(mixed_basis=False, body_to_body=False,
                         regular_wave=True)

    for name in ("float", "spar"):
        np.testing.assert_allclose(
            mixed.bodies[name].position, direct.bodies[name].position,
            rtol=0, atol=1e-10,
        )
        np.testing.assert_allclose(
            mixed.bodies[name].velocity, direct.bodies[name].velocity,
            rtol=0, atol=1e-10,
        )
    for field in ("stroke", "velocity", "force", "absorbed_power"):
        np.testing.assert_allclose(
            getattr(mixed.ptos["main"], field),
            getattr(direct.ptos["main"], field), rtol=0, atol=1e-7,
        )
    assert np.ptp(direct.wave_elevation) > .1
    assert np.max(np.abs(direct.bodies["spar"].position
                         - uncoupled.bodies["spar"].position)) > 1e-4

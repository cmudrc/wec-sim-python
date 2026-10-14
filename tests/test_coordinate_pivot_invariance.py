"""Equivalent public hinge coordinates must describe one physical motion."""

from pathlib import Path

import numpy as np

from wecsim import NoWave, WEC, WorldPoint


RM3 = Path(__file__).resolve().parents[1] / "examples/data/rm3.h5"
FLOAT_CENTER_Z = -0.72  # Pinned RM3 float CG; world hinge is at z=0.


def _run_pitch(*, local_pivot: bool, scale: int):
    wec = WEC("RM3 pitch coordinate conventions")
    body = wec.body("float", RM3, inertia=(0, 21_306_090.66, 0))
    pivot = body.at(0, 0, -FLOAT_CENTER_Z) if local_pivot else WorldPoint(0, 0, 0)
    pitch = wec.coordinate("pitch", body.move("pitch", scale=scale, pivot=pivot))
    wec.rotational_pto("hinge", pitch, stiffness=50_000, damping=10_000)
    return wec.run(
        NoWave(), dt=.01, end_time=.5, radiation_memory=.2,
        initial_coordinate={"pitch": .1 / scale},
    )


def test_pitch_pivot_and_coordinate_sign_preserve_motion_and_pto_work():
    reference = _run_pitch(local_pivot=False, scale=1)
    body = reference.bodies["float"]
    angle = body.position[:, 4]
    speed = body.velocity[:, 4]
    assert np.ptp(angle) > .02

    # Independent circular-hinge geometry, including the center's vertical
    # displacement and velocity that a constant small-angle map omits.
    np.testing.assert_allclose(body.position[:, 0],
                               FLOAT_CENTER_Z * np.sin(angle), atol=1e-12)
    np.testing.assert_allclose(body.position[:, 2],
                               FLOAT_CENTER_Z * np.cos(angle), atol=1e-12)
    np.testing.assert_allclose(body.velocity[:, 0],
                               FLOAT_CENTER_Z * np.cos(angle) * speed, atol=1e-12)
    np.testing.assert_allclose(body.velocity[:, 2],
                               -FLOAT_CENTER_Z * np.sin(angle) * speed, atol=1e-12)

    reference_pto = reference.ptos["hinge"]
    np.testing.assert_allclose(reference_pto.absorbed_power,
                               10_000 * reference_pto.velocity**2,
                               rtol=0, atol=1e-10)
    assert np.max(reference_pto.absorbed_power) > 10

    for local_pivot, scale in ((True, 1), (False, -1), (True, -1)):
        result = _run_pitch(local_pivot=local_pivot, scale=scale)
        for field in ("position", "velocity"):
            np.testing.assert_allclose(
                getattr(result.bodies["float"], field), getattr(body, field),
                rtol=0, atol=1e-11,
            )
        pto = result.ptos["hinge"]
        for field in ("stroke", "velocity", "force"):
            np.testing.assert_allclose(
                getattr(pto, field), getattr(reference_pto, field) / scale,
                rtol=0, atol=1e-11,
            )
        np.testing.assert_allclose(pto.absorbed_power,
                                   reference_pto.absorbed_power,
                                   rtol=0, atol=1e-10)

    # A nonunit coordinate scale changes PTO leverage, but the body center
    # must still follow the physical pitch angle on the exact circular path.
    scaled = _run_pitch(local_pivot=True, scale=2)
    scaled_body = scaled.bodies["float"]
    scaled_angle = scaled_body.position[:, 4]
    np.testing.assert_allclose(
        scaled_body.position[:, 0], FLOAT_CENTER_Z * np.sin(scaled_angle),
        rtol=0, atol=1e-12,
    )
    np.testing.assert_allclose(
        scaled_body.position[:, 2], FLOAT_CENTER_Z * np.cos(scaled_angle),
        rtol=0, atol=1e-12,
    )
    scaled_speed = scaled_body.velocity[:, 4]
    np.testing.assert_allclose(
        scaled_body.velocity[:, 0],
        FLOAT_CENTER_Z * np.cos(scaled_angle) * scaled_speed,
        rtol=0, atol=1e-12,
    )
    np.testing.assert_allclose(
        scaled_body.velocity[:, 2],
        -FLOAT_CENTER_Z * np.sin(scaled_angle) * scaled_speed,
        rtol=0, atol=1e-12,
    )
    np.testing.assert_allclose(2 * scaled.ptos["hinge"].stroke,
                               scaled_angle, rtol=0, atol=1e-12)

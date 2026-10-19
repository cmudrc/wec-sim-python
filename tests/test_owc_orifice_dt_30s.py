"""Keep the OWC source clock audit separate from published-case parity."""

from pathlib import Path

import numpy as np


FIXTURES = Path(__file__).parent / "fixtures/matlab"


def test_owc_source_surge_does_not_converge_monotonically_over_30s():
    with np.load(FIXTURES / "owc_orifice_dt_30s/source.npz",
                 allow_pickle=False) as saved, \
            np.load(FIXTURES / "owc_orifice_full_maxstep/source.npz",
                    allow_pickle=False) as full:
        time = saved["time"]
        interval = saved["interval"]
        motion = saved["motion"]
        velocity = saved["velocity"]
        wave = saved["wave"]
        components = saved["components"]
        work_internal = saved["work_internal"]
        work_output = saved["work_output"]

        assert time.shape == (6001,)
        np.testing.assert_allclose(time, np.arange(6001) * .005,
                                   rtol=0, atol=1e-11)
        np.testing.assert_array_equal(interval, [.005, .0025, .00125])
        assert motion.shape == velocity.shape == (3, 6001, 4)
        assert wave.shape == (3, 6001)
        assert components.shape == (3, 500, 4)
        assert np.isfinite(motion).all() and np.isfinite(velocity).all()

        # The first arm is the previously archived 130 s published-clock
        # source, not a new short-case baseline with a different startup.
        np.testing.assert_array_equal(motion[0], full["nominal_motion"][:6001])
        np.testing.assert_array_equal(wave[0], full["nominal_wave"][:6001])
        np.testing.assert_array_equal(components[0],
                                      full["nominal_components"])
        np.testing.assert_array_equal(
            wave[1:], np.broadcast_to(wave[0], wave[1:].shape))
        np.testing.assert_array_equal(
            components[1:],
            np.broadcast_to(components[0], components[1:].shape))

    first_surge_change = float(np.max(np.abs(motion[1, :, 0]
                                             - motion[0, :, 0])))
    second_surge_change = float(np.max(np.abs(motion[2, :, 0]
                                              - motion[1, :, 0])))
    # The smaller source step makes the surge self-change larger. This is a
    # pinned-source diagnostic and does not assert a Python parity failure.
    assert first_surge_change > .005
    assert second_surge_change > .03
    assert second_surge_change > 3 * first_surge_change

    assert np.all(work_internal > 0)
    assert np.all(np.abs(np.diff(work_internal)) / work_internal[:-1] < .002)
    # The internal orifice log has strictly increasing solver timestamps;
    # output-grid quadrature remains close enough to audit its work.
    assert np.all(np.abs(work_output - work_internal) / work_internal < .003)

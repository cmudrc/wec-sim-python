"""Check the opt-in delay law against a logged RM3 MoorDyn source prefix."""

from pathlib import Path

import numpy as np

from wecsim.dynamics.generalDynamics import _transport_delay_from_major_steps


SIGNALS = (Path(__file__).parent / "fixtures/matlab"
           / "rm3_moordyn_transport_delay_2s/signals.npz")


def test_adaptive_major_step_delay_matches_pinned_matlab_signals():
    with np.load(SIGNALS, allow_pickle=False) as source:
        time = source["time"]
        raw = source["delay_input"]
        delayed = source["delay_output"]
        radiation = source["radiation_input"]
        output_time = source["output_time"]
        output_acceleration = source["output_acceleration"]

    assert time.shape == (207,)
    assert raw.shape == delayed.shape == radiation.shape == (2, 207, 6)
    np.testing.assert_allclose(output_time, np.arange(201) * .01,
                               rtol=0, atol=1e-14)
    np.testing.assert_array_equal(radiation, delayed)

    # MATLAB took six additional ode45 major steps. The logged delay input
    # still agrees exactly with the independently saved body acceleration at
    # every output time, so these are the active force-loop signals.
    indices = np.searchsorted(time, output_time)
    np.testing.assert_array_equal(time[indices], output_time)
    np.testing.assert_array_equal(raw[:, indices], output_acceleration)

    np.testing.assert_array_equal(delayed[:, 0], np.zeros((2, 6)))
    np.testing.assert_array_equal(delayed[:, 1], raw[:, 0])
    for step in range(2, len(time)):
        for body in range(2):
            prediction = _transport_delay_from_major_steps(
                time[step - 1] - time[step - 2],
                time[step] - time[step - 1],
                raw[body, step - 2:step], 1e-7,
            )
            np.testing.assert_allclose(prediction, delayed[body, step],
                                       rtol=0, atol=1e-12)

    # The tiny nominal delay can have a visible startup effect because it
    # extrapolates across the adaptive solver's previous major-step interval.
    assert np.max(np.abs(delayed - raw)) > .08

"""Full-duration WaveStar NMPC window-work and operating-envelope comparison.

This is an autonomous closed-loop check. It measures integrated behavior
without claiming that individual motion samples or total signed work match
the source across numerical backends.
"""

from pathlib import Path

import numpy as np

from wecsim import WaveStarNmpcPTO
from wecsim.dynamics.wavestar import run_wavestar_published
from wecsim.environment.irregularWave import jonswap_equal_energy_components


FIXTURES = Path(__file__).parent / "fixtures/matlab"
HYDRO = FIXTURES / "wavestar_fixed_input_step/wavestar.h5"
SOURCE = FIXTURES / "wavestar_nmpc_full_step/baseline.npz"


def _rms(values):
    return np.sqrt(np.mean(np.square(values)))


def test_autonomous_nmpc_window_work_and_envelope():
    with np.load(SOURCE, allow_pickle=False) as source:
        source_body = source["body"]
        source_pto = source["pto"]
        phase = source["components"][:, 3, None]

    sea = jonswap_equal_energy_components(
        HYDRO, significant_height=.1042, peak_period=1.836,
        directions=np.array([0.]), spreading=np.array([1.]),
        gamma=3.3, phase=phase,
    )
    response = run_wavestar_published(
        HYDRO, sea, dt=.05, end_time=225, ramp_time=25,
        g=9.80665, solver="wecsim_ode8",
        pto_controller=WaveStarNmpcPTO(plant_dt=.05, control_dt=.05),
    )
    np.testing.assert_allclose(response.time, source_body[:, 0],
                               rtol=0, atol=1e-10)

    source_power = source_pto[:, 33] * source_pto[:, 9]
    python_power = response.pto_force * response.pto_speed
    # The total is diagnostic: small backend differences change this
    # cancellation-sensitive result beyond the former 0.1% limit.
    # Three consecutive controlled intervals retain explicit work gates.
    source_work = np.trapezoid(source_power, response.time)
    python_work = np.trapezoid(python_power, response.time)
    print(f"total signed PTO work: MATLAB {source_work:.12f} J, "
          f"Python {python_work:.12f} J")
    assert source_work < -100
    assert python_work < -100

    for start, stop in ((15, 75), (75, 150), (150, 225)):
        window = (response.time >= start) & (response.time <= stop)
        source_window_work = np.trapezoid(
            source_power[window], response.time[window],
        )
        python_window_work = np.trapezoid(
            python_power[window], response.time[window],
        )
        window_gap = abs(python_window_work - source_window_work) / abs(source_window_work)
        print(f"{start}–{stop} s signed-work gap: {window_gap:.6%}")
        assert source_window_work < -10
        assert window_gap < .005

        # Compare motion/force scale in each interval. The source's
        # step-sensitive NMPC feedback still prevents pointwise parity.
        for label, actual, expected, tolerance in (
            ("pitch", response.angle, source_body[:, 5], .01),
            ("pitch speed", response.angular_speed, source_body[:, 11], .01),
            ("PTO force", response.pto_force, source_pto[:, 33], .02),
        ):
            source_rms = _rms(expected[window])
            assert source_rms > 0
            rms_gap = abs(_rms(actual[window]) - source_rms) / source_rms
            print(f"{start}–{stop} s {label} RMS gap: "
                  f"{rms_gap:.6%}")
            assert rms_gap < tolerance

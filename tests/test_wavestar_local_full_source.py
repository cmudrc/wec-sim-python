"""Independent WaveStar NMPC motion against the local R2025b source run.

This is a published-input numerical comparison. The physical default and the
archived Linux source history are separate comparisons.
"""

from pathlib import Path

import numpy as np

from wecsim import WaveStarNmpcPTO
from wecsim.dynamics.wavestar import run_wavestar_published
from wecsim.environment.irregularWave import jonswap_equal_energy_components


FIXTURES = Path(__file__).parent / "fixtures/matlab"
HYDRO = FIXTURES / "wavestar_fixed_input_step/wavestar.h5"
LOCAL = FIXTURES / "wavestar_local_full_source/source.npz"
ARCHIVED = FIXTURES / "wavestar_nmpc_full_step/baseline.npz"


def test_full_225s_independent_ode8_against_local_r2025b():
    with np.load(LOCAL, allow_pickle=False) as source, \
            np.load(ARCHIVED, allow_pickle=False) as archived:
        time = source["time"]
        pitch = source["pitch"]
        speed = source["speed"]
        pto_force = source["pto_force"]
        pto_speed = source["pto_speed"]
        components = source["components"]
        wave = source["wave"]
        # The two MATLAB runs use identical hydrodynamics and incident sea;
        # their later closed-loop trajectories are platform-sensitive.
        np.testing.assert_array_equal(components, archived["components"])
        np.testing.assert_array_equal(wave, archived["wave"])
        assert np.max(np.abs(pitch - archived["body"][:, 5])) > .007

    assert all(values.shape == (4501,) for values in (
        time, pitch, speed, pto_force, pto_speed,
    ))
    np.testing.assert_allclose(time, np.arange(4501) * .05, rtol=0, atol=1e-10)
    assert all(np.isfinite(values).all() for values in (
        time, pitch, speed, pto_force, pto_speed,
    ))
    sea = jonswap_equal_energy_components(
        HYDRO, significant_height=.1042, peak_period=1.836,
        directions=np.array([0.]), spreading=np.array([1.]),
        gamma=3.3, phase=components[:, 3, None],
    )
    response = run_wavestar_published(
        HYDRO, sea, dt=.05, end_time=225, ramp_time=25,
        g=9.80665, solver="wecsim_ode8",
        pto_controller=WaveStarNmpcPTO(plant_dt=.05, control_dt=.05),
    )
    np.testing.assert_allclose(response.time, time, rtol=0, atol=1e-10)
    assert np.max(np.abs(response.angle - pitch)) < .012
    assert np.max(np.abs(response.angular_speed - speed)) < .08
    assert np.max(np.abs(response.pto_force - pto_force)) < 5.
    python_work = np.trapezoid(
        response.pto_force * response.pto_speed, time,
    )
    source_work = np.trapezoid(pto_force * pto_speed, time)
    assert abs(python_work - source_work) / abs(source_work) < .002

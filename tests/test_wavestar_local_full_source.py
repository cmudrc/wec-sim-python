"""Independent WaveStar NMPC motion against the local R2025b source run.

This is a published-input numerical comparison. The physical default and the
archived Linux source history are separate comparisons.
"""

from pathlib import Path
import platform

import numpy as np
import pytest

from wecsim import WaveStarNmpcPTO
from wecsim.dynamics.wavestar import run_wavestar_published
from wecsim.environment.irregularWave import jonswap_equal_energy_components


FIXTURES = Path(__file__).parent / "fixtures/matlab"
HYDRO = FIXTURES / "wavestar_fixed_input_step/wavestar.h5"
LOCAL = FIXTURES / "wavestar_local_full_source/source.npz"
ARCHIVED = FIXTURES / "wavestar_nmpc_full_step/baseline.npz"


@pytest.mark.skipif(
    platform.system() != "Darwin" or platform.machine() != "arm64",
    reason="the full-run source baseline was generated on macOS ARM64",
)
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
    pitch_gap = np.abs(response.angle - pitch)
    speed_gap = np.abs(response.angular_speed - speed)
    force_gap = np.abs(response.pto_force - pto_force)
    python_work = np.trapezoid(
        response.pto_force * response.pto_speed, time,
    )
    source_work = np.trapezoid(pto_force * pto_speed, time)
    work_gap = abs(python_work - source_work) / abs(source_work)
    print("WaveStar numerical runner:", platform.platform(), platform.machine())
    np.show_runtime()
    print("Full WaveStar source gaps:", float(pitch_gap.max()),
          float(speed_gap.max()), float(force_gap.max()), float(work_gap))
    for end in (15, 20, 30, 75, 150, 225):
        count = round(end / .05) + 1
        print("Pitch prefix", end, float(pitch_gap[:count].max()))
    assert pitch_gap.max() < .012
    assert speed_gap.max() < .08
    assert force_gap.max() < 5.
    assert work_gap < .002

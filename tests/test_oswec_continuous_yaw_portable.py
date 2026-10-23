"""Full-duration physical parity for the repository-owned OSWEC yaw control."""

import hashlib
import json
from pathlib import Path

import numpy as np
import pytest
from scipy.io import loadmat

from tools.oswec_continuous_yaw_step_audit import audit
from wecsim import PMWave, WEC, WorldPoint


FIXTURES = Path(__file__).resolve().parent / "fixtures/matlab"
SOURCE = FIXTURES / "oswec_continuous_yaw"
HYDRO = FIXTURES / "oswec_multi_wave/oswec.h5"


@pytest.mark.parametrize("seed", (1, 2, 3))
def test_adjusted_250s_source_is_step_stable(seed):
    """The archived corrected source must satisfy the predeclared limits."""
    manifest = json.loads((SOURCE / "manifest.json").read_text())
    patch = FIXTURES / "oswec_held_yaw/continuous_heading.patch"
    assert (manifest["continuous_heading_patch_sha256"]
            == hashlib.sha256(patch.read_bytes()).hexdigest())
    report = audit(SOURCE, seed)
    assert report["within_limits"]
    assert report == manifest["seeds"][str(seed)]


@pytest.mark.parametrize("seed", (1, 2, 3))
def test_adjusted_250s_public_physical_path(seed):
    """Advance each sea from phases and HDF5, without source states or forces."""
    source = loadmat(SOURCE / f"seed-{seed}-dt-0.005.mat", squeeze_me=True)
    wec = WEC("OSWEC continuous-heading physical parity")
    flap = wec.body(
        "flap", HYDRO, mass=12700, inertia=(1.85e6,) * 3,
        passive_yaw=True, passive_yaw_threshold=0,
    )
    wec.body("base", HYDRO, mass=999, inertia=(999,) * 3)
    yaw = wec.coordinate("yaw", flap.move("yaw", pivot=WorldPoint(0, 0, -8.9)))
    wec.rotational_pto("hinge", yaw, damping=120000)
    result = wec.run(
        PMWave(2.5, 8, direction=10, seed=seed, phase_generator="matlab"),
        dt=0.005, end_time=250, ramp_time=100, radiation_memory=40,
    )
    time = np.asarray(source["time"]).ravel()
    source_yaw = np.asarray(source["yaw"]).ravel()
    source_speed = np.asarray(source["yawRate"]).ravel()
    source_torque = np.asarray(source["ptoForce"]).ravel()
    source_power = np.asarray(source["ptoPower"]).ravel()
    source_wave = np.asarray(source["wave"])[:, 1]
    np.testing.assert_allclose(result.time, time, rtol=0, atol=1e-10)
    assert np.max(np.abs(result.wave_elevation - source_wave)) < 1e-11
    assert np.max(np.abs(result.coordinates["yaw"].position - source_yaw)) < 1e-4
    assert np.max(np.abs(result.coordinates["yaw"].velocity - source_speed)) < 2e-5
    assert np.max(np.abs(result.ptos["hinge"].force - source_torque)) < 2.4
    source_work = -np.trapezoid(source_power, time)
    physical_work = np.trapezoid(result.ptos["hinge"].absorbed_power, time)
    assert source_work > 0
    assert abs(physical_work / source_work - 1) < 0.0002

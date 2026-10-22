"""Full-duration physical parity against a source-stable continuous-yaw step."""

import os
from pathlib import Path

import numpy as np
import pytest
from scipy.io import loadmat

from wecsim import PMWave, WEC, WorldPoint


APPLICATIONS = os.environ.get("WEC_SIM_APPLICATIONS_DIR")
REFERENCE = os.environ.get("WEC_SIM_CONTINUOUS_REFINEMENT_DIR")
PHASE_SEED = os.environ.get("WEC_SIM_PHASE_SEED")
pytestmark = pytest.mark.skipif(
    not (APPLICATIONS and REFERENCE and PHASE_SEED),
    reason="pinned continuous-yaw source refinement not provided",
)


def test_continuous_heading_physical_full_duration():
    seed = int(PHASE_SEED)
    assert seed in (1, 2, 3)
    hydro = (Path(APPLICATIONS)
             / "_Common_Input_Files/OSWEC/hydroData/oswec.h5")
    source = loadmat(
        Path(REFERENCE) / f"seed-{seed}-dt-0.005.mat", squeeze_me=True,
    )
    wec = WEC("OSWEC continuous-heading physical parity")
    flap = wec.body(
        "flap", hydro, mass=12700, inertia=(1.85e6,) * 3,
        passive_yaw=True, passive_yaw_threshold=0,
    )
    wec.body("base", hydro, mass=999, inertia=(999,) * 3)
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
    print(f"seed {seed} continuous 250 s: source work {source_work:.6f} J, "
          f"physical work {physical_work:.6f} J")

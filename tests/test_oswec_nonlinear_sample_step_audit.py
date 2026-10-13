"""Separate OSWEC nonlinear-force sampling from fixed-step ode4 integration."""

import os
from pathlib import Path

import numpy as np
import pytest
from scipy.io import loadmat

from wecsim import RegularWave, WEC, WorldPoint


FIXTURE = (Path(__file__).parent / "fixtures" / "matlab"
           / "oswec_nonlinear_viz_published_motion.npz")


@pytest.mark.skipif(
    not (os.getenv("WEC_SIM_OSWEC_SAMPLE_STEP_DIR")
         and os.getenv("WEC_SIM_OSWEC_H5")
         and os.getenv("WEC_SIM_OSWEC_FLAP_STL")),
    reason="pinned OSWEC nonlinear sample-step audit assets not provided",
)
def test_nonlinear_force_sample_step_at_fixed_solver_step():
    directory = Path(os.environ["WEC_SIM_OSWEC_SAMPLE_STEP_DIR"])
    fine = loadmat(directory / "fine.mat")
    coarse_force = loadmat(directory / "coarse-force.mat")
    with np.load(FIXTURE) as source:
        published = {key: source[key] for key in source.files}

    for arm, expected_nonlinear_dt in ((fine, .0125), (coarse_force, .1)):
        assert float(arm["nonlinearDt"].item()) == expected_nonlinear_dt
        np.testing.assert_allclose(arm["time"].ravel()[::8],
                                   published["time"].ravel(),
                                   rtol=0, atol=1e-12)
        for key in ("pose", "velocity", "acceleration", "forceExcitation",
                    "forceRadiationDamping", "forceAddedMass", "forceRestoring"):
            assert arm[key].shape == (9601, 6)
            assert np.isfinite(arm[key]).all()

    # This is the previously measured source refinement, held as an identity
    # check before interpreting the new force-sampling arm. The split arm is
    # diagnostic: either outcome changes the next physical investigation.
    fine_pose = fine["pose"][::8]
    split_pose = coarse_force["pose"][::8]
    coarse_pose = published["pose"]
    fine_pitch_gap = np.max(np.abs(fine_pose[:, 4] - coarse_pose[:, 4]))
    assert 0.0070 < fine_pitch_gap < 0.0074

    for label, key, column in (
        ("pitch", "pose", 4), ("surge", "pose", 0),
        ("heave", "pose", 2), ("pitch speed", "velocity", 4),
        ("restoring pitch moment", "forceRestoring", 4),
        ("excitation pitch moment", "forceExcitation", 4),
        ("added-mass pitch moment", "forceAddedMass", 4),
    ):
        coarse = published[key][:, column]
        refined = fine[key][::8, column]
        split = coarse_force[key][::8, column]
        published_to_fine = np.max(np.abs(refined - coarse))
        published_to_split = np.max(np.abs(split - coarse))
        split_to_fine = np.max(np.abs(refined - split))
        denominator = np.sum((refined - coarse) ** 2)
        captured = (np.sum((split - coarse) * (refined - coarse))
                    / denominator if denominator else np.nan)
        print(f"{label}: published-to-fine={published_to_fine:.9g}, "
              f"published-to-split={published_to_split:.9g}, "
              f"split-to-fine={split_to_fine:.9g}, "
              f"solver-only projection={captured:.6g}")

    # The public physical solver is an independent control against the same
    # fine source arm. It uses only the HDF5, STL, regular sea, and case inputs.
    wec = WEC("OSWEC nonlinear sample-step audit")
    flap = wec.body(
        "flap", os.environ["WEC_SIM_OSWEC_H5"], mass=127000,
        inertia=(1.85e6,) * 3,
        geometry_file=os.environ["WEC_SIM_OSWEC_FLAP_STL"],
        nonlinear_hydro="instantaneous",
    )
    base = wec.fixed_body(
        "base", center_gravity=(0, 0, -10.9),
        mass=999, inertia=(1, 1, 1),
    )
    wec.fixed_hinge(
        flap, base, location=WorldPoint(0, 0, -10),
        pto_location=WorldPoint(0, 0, -8.9),
    )
    response = wec.run(RegularWave(2.5, 8), dt=.1, end_time=120,
                       ramp_time=40)
    pose = response.bodies["flap"].position
    speed = response.bodies["flap"].velocity
    for axis, limit in ((4, 2e-4), (0, 1e-3), (2, 4e-4)):
        assert np.max(np.abs(pose[:, axis] - fine_pose[:, axis])) < limit
    for axis, limit in ((4, 1.5e-4), (0, 8e-4), (2, 3.5e-4)):
        assert np.max(np.abs(speed[:, axis]
                             - fine["velocity"][::8, axis])) < limit
    print("Python pitch gap to source with coarse nonlinear sampling:",
          np.max(np.abs(pose[:, 4] - split_pose[:, 4])))

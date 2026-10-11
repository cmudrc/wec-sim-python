"""Independent nonlinear OSWEC flap motion against refined MATLAB runs."""

import os
from pathlib import Path

import numpy as np
import pytest
from scipy.io import loadmat

from wecsim import RegularWave, WEC, WorldPoint


@pytest.mark.skipif(
    not (os.getenv("WEC_SIM_OSWEC_REFINEMENT_DIR")
         and os.getenv("WEC_SIM_OSWEC_H5")
         and os.getenv("WEC_SIM_OSWEC_FLAP_STL")),
    reason="pinned refined MATLAB motion, HDF5, and flap STL not provided",
)
def test_published_nonlinear_oswec_hinge_against_refined_matlab():
    """Run Python from case inputs, without replaying any MATLAB body states."""
    reference = Path(os.environ["WEC_SIM_OSWEC_REFINEMENT_DIR"])
    coarse = loadmat(reference / "dt-0.1.mat")
    medium = loadmat(reference / "dt-0.05.mat")
    fine = loadmat(reference / "dt-0.025.mat")
    finest = loadmat(reference / "dt-0.0125.mat")
    for source, dt in ((coarse, 0.1), (medium, 0.05),
                       (fine, 0.025), (finest, 0.0125)):
        assert float(source["dt"].item()) == dt
        assert float(source["nonlinearDt"].item()) == dt
        assert str(source["solver"].item()) == "ode4"
        assert len(source["time"].ravel()) == round(120 / dt) + 1
    np.testing.assert_allclose(coarse["time"].ravel(),
                               medium["time"].ravel()[::2], rtol=0, atol=1e-12)
    np.testing.assert_allclose(medium["time"].ravel(),
                               fine["time"].ravel()[::2], rtol=0, atol=1e-12)
    np.testing.assert_allclose(fine["time"].ravel(),
                               finest["time"].ravel()[::2], rtol=0, atol=1e-12)
    coarse_change = np.max(np.abs(
        coarse["pose"][:, 4] - medium["pose"][::2, 4]
    ))
    refined_change = np.max(np.abs(
        medium["pose"][:, 4] - fine["pose"][::2, 4]
    ))
    finest_change = np.max(np.abs(
        fine["pose"][:, 4] - finest["pose"][::2, 4]
    ))
    assert 0.005 < coarse_change < 0.006
    assert refined_change < 0.3 * coarse_change
    assert finest_change < 0.3 * refined_change

    wec = WEC("Published nonlinear OSWEC")
    flap = wec.body(
        "flap", os.environ["WEC_SIM_OSWEC_H5"], mass=127000,
        inertia=(1.85e6, 1.85e6, 1.85e6),
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
    result = wec.run(
        RegularWave(2.5, 8), dt=0.1, end_time=120, ramp_time=40,
    )
    np.testing.assert_allclose(result.time, finest["time"].ravel()[::8],
                               rtol=0, atol=1e-12)
    position = result.bodies["flap"].position
    velocity = result.bodies["flap"].velocity
    source_position = finest["pose"][::8]
    source_velocity = finest["velocity"][::8]
    assert np.max(np.abs(position[:, 4] - source_position[:, 4])) < 2e-4
    assert np.max(np.abs(velocity[:, 4] - source_velocity[:, 4])) < 1.5e-4
    assert np.max(np.abs(position[:, 0] - source_position[:, 0])) < 1e-3
    assert np.max(np.abs(position[:, 2] - source_position[:, 2])) < 4e-4
    assert np.max(np.abs(velocity[:, 0] - source_velocity[:, 0])) < 8e-4
    assert np.max(np.abs(velocity[:, 2] - source_velocity[:, 2])) < 3.5e-4
    np.testing.assert_allclose(result.ptos["hinge"].force, 0, rtol=0, atol=0)

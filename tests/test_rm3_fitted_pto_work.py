"""Numerical PTO-work gate for the pinned nonpassive RM3 radiation fit."""

from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import trapezoid

from wecsim import RegularCICWave, WEC


ROOT = Path(__file__).resolve().parents[1]
HYDRO = ROOT / "examples/data/rm3.h5"
SOURCE = ROOT / "tests/fixtures/matlab/rm3_b2b_cic/source.npz"


@pytest.mark.parametrize("case,body_to_body", [
    ("B2B_Case5", False), ("B2B_Case6", True),
])
def test_opt_in_source_fit_pto_work(case, body_to_body):
    """Pair source numerics without changing the physical radiation default."""
    with np.load(SOURCE, allow_pickle=False) as saved:
        source_pto = saved[f"{case}_pto1"].copy()
    assert source_pto.shape == (4001, 25)

    wec = WEC(case, body_to_body=body_to_body)
    float_body = wec.body("float", HYDRO, inertia=(0, 21_306_090.66, 0))
    spar = wec.body("spar", HYDRO, inertia=(0, 94_407_091.24, 0))
    wec.floating_joint(float_body, spar, damping=1_200_000,
                       radiation_method="state_space")
    result = wec.run(RegularCICWave(2.5, 8), dt=.1, end_time=400,
                     ramp_time=100, radiation_memory=60)
    assert result.case["simulation"]["radiation_method"] == "state_space"
    np.testing.assert_allclose(result.time, source_pto[:, 0], rtol=0,
                               atol=1e-10)

    pto = result.ptos["relative_heave"]
    assert np.max(np.abs(pto.absorbed_power + source_pto[:, 21])) < 6_000
    source_work = -trapezoid(source_pto[:, 21], x=source_pto[:, 0])
    python_work = trapezoid(pto.absorbed_power, x=result.time)
    assert source_work > 0
    assert abs(python_work - source_work) / source_work < .01

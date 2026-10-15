"""Full-duration public WEC gates for the RM3 convolution B2B cases."""

from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import trapezoid

from wecsim import RegularCICWave, WEC


ROOT = Path(__file__).resolve().parents[1]
HYDRO = ROOT / "examples/data/rm3.h5"
SOURCE = ROOT / "tests/fixtures/matlab/rm3_b2b_cic/source.npz"


@pytest.mark.parametrize("case,body_to_body", [
    ("B2B_Case3", False), ("B2B_Case4", True),
])
def test_public_rm3_convolution_motion_and_pto_work(case, body_to_body):
    """Advance the published input without source motion or force histories."""
    with np.load(SOURCE, allow_pickle=False) as saved:
        body1 = saved[f"{case}_body1"].copy()
        body2 = saved[f"{case}_body2"].copy()
        source_pto = saved[f"{case}_pto1"].copy()
    assert body1.shape == body2.shape == source_pto.shape == (4001, 25)

    wec = WEC(case, body_to_body=body_to_body)
    float_body = wec.body("float", HYDRO, inertia=(0, 21_306_090.66, 0))
    spar = wec.body("spar", HYDRO, inertia=(0, 94_407_091.24, 0))
    wec.floating_joint(float_body, spar, damping=1_200_000)
    response = wec.run(
        RegularCICWave(2.5, 8), dt=.1, end_time=400,
        ramp_time=100, radiation_memory=60,
    )
    assert response.case.get("body_to_body", False) is body_to_body
    assert response.case["wave"]["type"] == "regularCIC"
    assert response.case["pto"]["damping"] == 1_200_000
    np.testing.assert_allclose(response.time, source_pto[:, 0],
                               rtol=0, atol=1e-10)

    for name, expected in (("float", body1), ("spar", body2)):
        body = response.bodies[name]
        for dof, position_limit, speed_limit in (
            (0, .075, .04), (2, .005, .0035), (4, .002, .0011),
        ):
            assert np.max(np.abs(body.position[:, dof]
                                 - expected[:, 1 + dof])) < position_limit
            assert np.max(np.abs(body.velocity[:, dof]
                                 - expected[:, 7 + dof])) < speed_limit

    pto = response.ptos["relative_heave"]
    assert np.max(np.abs(pto.stroke - source_pto[:, 3])) < .008
    assert np.max(np.abs(pto.velocity - source_pto[:, 9])) < .005
    assert np.max(np.abs(pto.force - source_pto[:, 15])) < 6_000
    assert np.max(np.abs(pto.absorbed_power + source_pto[:, 21])) < 7_000
    source_work = -trapezoid(source_pto[:, 21], x=source_pto[:, 0])
    python_work = trapezoid(pto.absorbed_power, x=response.time)
    assert source_work > 0
    assert abs(python_work - source_work) / source_work < .01

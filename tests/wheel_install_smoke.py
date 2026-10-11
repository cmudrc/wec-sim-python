"""Run a configured two-body WEC using only an installed wheel and HDF5 input."""

from pathlib import Path
import sys

import numpy as np

import wecsim
from wecsim import NoWave, WEC


hydro = Path(sys.argv[1])
assert Path(wecsim.__file__).is_relative_to(Path(sys.prefix))

wec = WEC("Wheel install smoke")
float_body = wec.body("float", hydro)
spar_body = wec.body("spar", hydro)
wec.coordinate("float_heave", float_body.move("heave"))
wec.coordinate("spar_heave", spar_body.move("heave"))
wec.pto(
    "main", float_body.at(2, 0, 0.72), spar_body.at(-1, 0, 21.29),
    axis=(0, 0, 1), stiffness=1_000, damping=1_200,
)

result = wec.run(
    NoWave(), dt=0.1, end_time=0.2, radiation_memory=0.2,
    initial_coordinate={"float_heave": 0.1},
)
pto = result.ptos["main"]
np.testing.assert_array_equal(result.time, [0, 0.1, 0.2])
np.testing.assert_allclose(pto.stroke[0], -0.1, rtol=0, atol=1e-12)
np.testing.assert_allclose(
    pto.force, -1_000 * pto.stroke - 1_200 * pto.velocity,
    rtol=0, atol=1e-9,
)
np.testing.assert_allclose(
    pto.absorbed_power, 1_200 * pto.velocity**2,
    rtol=0, atol=1e-9,
)
assert np.isfinite(result.bodies["float"].position).all()
assert np.isfinite(result.bodies["spar"].position).all()

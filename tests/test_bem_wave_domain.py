"""Public wave inputs must stay inside their hydrodynamic coefficient grid."""

from pathlib import Path

import numpy as np
import pytest

from wecsim import PMWave, RegularWave, WEC


SPHERE = Path(__file__).resolve().parents[1] / "examples/data/sphere.h5"
OSWEC = Path(__file__).resolve().parents[1] / "examples/data/oswec.h5"


def _heave_wec(hydro_file):
    wec = WEC("wave coverage")
    body = wec.body("body", hydro_file)
    wec.coordinate("heave", body.move("heave"))
    return wec


@pytest.mark.parametrize("wave, settings, error", (
    (RegularWave(.5, 8, direction=15),
     {"ramp_time": 0}, "wave direction"),
    (PMWave(.5, 8, direction=15, seed=1, frequency_count=16),
     {"radiation_memory": .2}, "wave direction"),
    (RegularWave(.5, .5),
     {"ramp_time": 0}, "wave frequency"),
))
def test_single_heading_hdf5_rejects_unsupported_incidence(wave, settings, error):
    """MATLAB checks the BEM range during setup; Python must not return zero force."""
    wec = _heave_wec(SPHERE)
    with pytest.raises(ValueError, match=error):
        wec.run(wave, dt=.1, end_time=.2, **settings)

    # A rejected run leaves the same Python configuration usable with a
    # supported incident wave and a nonzero hydrodynamic response.
    result = wec.run(RegularWave(.5, 8), dt=.1, end_time=2, ramp_time=0)
    assert np.ptp(result.bodies["body"].position[:, 2]) > .01


def test_multidirectional_hdf5_rejects_heading_extrapolation():
    wec = _heave_wec(OSWEC)
    with pytest.raises(ValueError, match="wave direction"):
        wec.run(RegularWave(.5, 8, direction=95),
                dt=.1, end_time=.2, ramp_time=0)
    result = wec.run(RegularWave(.5, 8, direction=15),
                     dt=.1, end_time=.2, ramp_time=0)
    assert np.isfinite(result.bodies["body"].position).all()

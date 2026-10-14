"""A configured barge mass must change its still-water heave response."""

from pathlib import Path

import numpy as np

from wecsim import RegularWave, WEC


HYDRO = Path(__file__).parent / "fixtures/matlab/gbm_barge.h5"
SOURCE = Path(__file__).parent / "fixtures/matlab/gbm_custom_mass/heavy_heave.csv"


def _still_water(mass):
    wec = WEC("GBM barge static mass")
    body = wec.body("barge", HYDRO, mass=mass,
                    inertia=(6.667e7, 2.167e9, 2.167e9))
    wec.floating_gbm(body)
    return wec.run(RegularWave(height=0, period=8), dt=0.05,
                   end_time=1, ramp_time=0)


def test_custom_mass_heave_against_local_matlab():
    source = np.loadtxt(SOURCE, delimiter=",")
    assert source.shape == (21, 3)
    result = _still_water(4_100_000)
    assert np.allclose(result.time, source[:, 0], rtol=0, atol=1e-12)
    heave = result.bodies["barge"].position[:, 2]
    speed = result.bodies["barge"].velocity[:, 2]
    assert np.isfinite(heave).all() and np.isfinite(speed).all()
    assert np.max(np.abs(heave - source[:, 1])) < 0.003
    assert np.max(np.abs(speed - source[:, 2])) < 0.009
    assert heave[-1] < -0.04


def test_equilibrium_mass_has_no_still_water_motion():
    result = _still_water("equilibrium")
    assert np.all(result.bodies["barge"].position == 0)
    assert np.all(result.bodies["barge"].velocity == 0)

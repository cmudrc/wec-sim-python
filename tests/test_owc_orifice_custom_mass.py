"""Pair the OWC's still-water mass offset with a bounded MATLAB run."""

from pathlib import Path

import numpy as np

from wecsim import OrificePTO, PMWave, WEC


FIXTURES = Path(__file__).parent / "fixtures/matlab"
HYDRO = FIXTURES / "owc_orifice.h5"
SOURCE = FIXTURES / "owc_orifice_custom_mass/heave.csv"


def _run(mass):
    with np.load(FIXTURES / "owc_orifice_prefix.npz", allow_pickle=False) as saved:
        orifice = OrificePTO(*saved["parameters"])
    wec = WEC("OWC custom mass")
    body = wec.body("OWC", HYDRO, mass=mass,
                    inertia=(99.28, 11.04, 99.2))
    wec.floating_gbm(
        body, orifice=orifice, orifice_force_path="published_owc",
        heave_linear_damping=100, mode_linear_damping=100,
        heave_drag_cd=1.2, heave_drag_area=8,
        pitch_drag_cd=1.2, pitch_drag_area=8,
    )
    return wec.run(
        PMWave(height=1e-8, period=4, seed=128, phase_generator="matlab"),
        dt=.005, end_time=.05, ramp_time=10, radiation_memory=15,
    )


def test_custom_mass_heave_against_local_matlab():
    source = np.loadtxt(SOURCE, delimiter=",")
    assert source.shape == (11, 3)
    result = _run(689.049)
    np.testing.assert_allclose(result.time, source[:, 0], rtol=0, atol=1e-12)
    heave = result.bodies["OWC"].position[:, 2]
    speed = result.bodies["OWC"].velocity[:, 2]
    assert np.isfinite(heave).all() and np.isfinite(speed).all()
    assert np.max(np.abs(heave - source[:, 1])) < 5e-6
    assert np.max(np.abs(speed - source[:, 2])) < 6e-4
    assert heave[-1] < -5e-4


def test_equilibrium_mass_has_negligible_heave():
    source = np.loadtxt(SOURCE.with_name("equilibrium_heave.csv"), delimiter=",")
    assert source.shape == (11, 3)
    result = _run("equilibrium")
    np.testing.assert_allclose(result.time, source[:, 0], rtol=0, atol=1e-12)
    np.testing.assert_allclose(result.bodies["OWC"].position[:, 2],
                               source[:, 1], rtol=0, atol=1e-9)
    np.testing.assert_allclose(result.bodies["OWC"].velocity[:, 2],
                               source[:, 2], rtol=0, atol=1e-8)

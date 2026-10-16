"""Public wave inputs must stay inside their hydrodynamic coefficient grid."""

from pathlib import Path

import h5py
import numpy as np
import pytest

from wecsim import PMWave, RegularWave, WEC
from wecsim.environment.irregularWave import (
    FullDirectionalComponents, synthesize_full_directional_response,
)


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


@pytest.mark.parametrize("interpolation", ("linear", "spline_frequency"))
@pytest.mark.parametrize("edge", ("below", "above"))
def test_full_directional_force_rejects_frequency_extrapolation(interpolation, edge):
    """A supplied spectrum cannot bypass MATLAB's BEM-frequency domain check."""
    with h5py.File(OSWEC) as hydro:
        bem_omega = np.asarray(hydro["simulation_parameters/w"]).ravel()
    components = FullDirectionalComponents(
        omega=(np.array([bem_omega[0] - .1, bem_omega[100]]) if edge == "below"
               else np.array([bem_omega[100], bem_omega[-1] + .1])),
        spectral_amplitude=np.full((2, 2), .01),
        d_omega=np.full(2, .1), directions=np.array([0., 10.]),
        d_theta=np.full(2, .1), phase=np.zeros((2, 2)),
    )
    settings = dict(dt=.1, end_time=.1, ramp_time=0,
                    excitation_interpolation=interpolation)
    with pytest.raises(ValueError, match="full-directional frequencies"):
        synthesize_full_directional_response(OSWEC, components, **settings)

    supported = FullDirectionalComponents(
        omega=np.array([bem_omega[100], bem_omega[200]]),
        spectral_amplitude=components.spectral_amplitude,
        d_omega=components.d_omega, directions=components.directions,
        d_theta=components.d_theta, phase=components.phase,
    )
    response = synthesize_full_directional_response(OSWEC, supported, **settings)
    assert np.isfinite(response.excitation_force).all()
    assert np.max(np.abs(response.excitation_force)) > 0

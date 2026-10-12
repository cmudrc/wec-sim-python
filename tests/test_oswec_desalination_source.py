"""Published Desalination sea and linkage, without claiming fluid parity."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim import PitchRodLinkage
from wecsim.environment.irregularWave import (
    pm_equal_energy_components, synthesize_irregular_response,
)
from wecsim.hydrodynamics.morison import MorisonElement, irregular_morison_source_drag


REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
OSWEC_H5 = os.environ.get("WEC_SIM_OSWEC_H5")
pytestmark = pytest.mark.skipif(
    not (REFERENCE and OSWEC_H5),
    reason="pinned OSWEC Desalination MATLAB output and HDF5 absent",
)


def _source(name):
    values = np.loadtxt(
        Path(REFERENCE) / f"OSWEC_DESALINATION_SOURCE_{name}.csv",
        delimiter=",",
    )
    assert values.ndim == 2 and np.isfinite(values).all(), name
    return values


def test_published_250_bin_wave():
    source = _source("components")
    assert source.shape == (250, 4)
    components = pm_equal_energy_components(
        OSWEC_H5, significant_height=2.64, peak_period=9.86,
        directions=(0,), spreading=(1,), count=250,
        phase=source[:, 3:4],
    )
    actual = np.column_stack((
        components.omega, components.spectral_amplitude,
        components.d_omega, components.phase.ravel(),
    ))
    np.testing.assert_allclose(actual, source, rtol=0, atol=1e-12)
    sea = synthesize_irregular_response(
        OSWEC_H5, components, dt=.01, end_time=300, ramp_time=50,
    )
    wave = _source("wave")
    assert wave.shape == (30_001, 2)
    np.testing.assert_allclose(sea.elevation, wave[:, 1], rtol=0,
                               atol=2e-12)
    body = _source("Desalination_body1")
    np.testing.assert_allclose(sea.excitation_force[:, :6], body[:, 19:25],
                               rtol=0, atol=1e-6)


def test_published_body_local_rod_motion():
    body = _source("Desalination_body1")
    base = _source("Desalination_body2")
    pto = _source("Desalination_pto1")
    measured = _source("simout")
    assert all(x.shape[0] == 30_001 for x in (body, base, pto, measured))
    angle, angular_speed = body[:, 5], body[:, 11]
    np.testing.assert_allclose(body[:, 1], 5 * np.sin(angle),
                               rtol=0, atol=1e-9)
    np.testing.assert_allclose(body[:, 3],
                               -3.9 + 5 * (np.cos(angle) - 1),
                               rtol=0, atol=1e-9)
    assert np.max(np.abs(base[:, 1:7] - base[0, 1:7])) < 1e-12
    linkage = PitchRodLinkage(
        anchor=(5.6021271782, -8.7), hinge=(0, -8.9),
        body_center=(0, -3.9), body_point=(.9, -3.1),
    )
    stroke, jacobian = linkage.stroke_and_jacobian(angle)
    np.testing.assert_allclose(stroke, pto[:, 3], rtol=0, atol=1e-9)
    np.testing.assert_allclose(jacobian * angular_speed, pto[:, 9],
                               rtol=0, atol=1e-10)
    np.testing.assert_allclose(measured[:, 1], pto[:, 9],
                               rtol=0, atol=1e-12)


def test_published_hydraulic_measurements_are_active():
    mechanical = _source("simout")
    hydraulic = _source("simout1")
    pto = _source("Desalination_pto1")
    assert mechanical.shape == (30_001, 4)
    assert hydraulic.shape == (30_001, 8)
    assert pto.shape == (30_001, 49)
    np.testing.assert_allclose(mechanical[:, 3],
                               mechanical[:, 1] * mechanical[:, 2],
                               rtol=0, atol=1e-7)
    # The source PTO actuation block delays the measured cylinder force
    # by one discrete output interval to break its algebraic loop.
    np.testing.assert_allclose(pto[1:, 39], mechanical[:-1, 2],
                               rtol=0, atol=1e-8)
    assert abs(pto[0, 39]) < 1e-12
    assert np.max(np.abs(mechanical[:, 2])) > 1e6
    assert np.max(hydraulic[:, 6]) > 5e6
    assert np.max(hydraulic[:, 2:6]) > .1


def test_published_flap_morison_drag_on_saved_body_states():
    body = _source("Desalination_body1")
    source = _source("flap_forces")
    wave = _source("components")
    components = pm_equal_energy_components(
        OSWEC_H5, significant_height=2.64, peak_period=9.86,
        directions=(0,), spreading=(1,), count=250,
        phase=wave[:, 3:4],
    )
    elements = [MorisonElement(
        point=(0, 0, z), drag_coefficient=(1, 1, 1),
        added_mass_coefficient=(0, 0, 0), area=(32.4, 0, 32.4),
        volume=0,
    ) for z in (-3, -1.2, .6, 2.4, 4.2)]
    samples = np.unique(np.r_[
        np.linspace(0, 30_000, 601, dtype=int),
        np.argmax(np.abs(source[:, 1:7]), axis=0),
    ])
    actual = np.array([irregular_morison_source_drag(
        elements, time=body[i, 0], position=body[i, 1:7],
        velocity=body[i, 7:13], components=components,
        water_depth=10.9, ramp_time=50,
    ) for i in samples])
    # WEC-Sim logs this body channel with the opposite sign of the source
    # Morison function's applied force.
    np.testing.assert_allclose(actual, -source[samples, 1:7],
                               rtol=0, atol=1e-6)

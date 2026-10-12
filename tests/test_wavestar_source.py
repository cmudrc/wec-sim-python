"""Published WaveStar linkage and unforced PTO against pinned MATLAB output."""

import os
from pathlib import Path

import numpy as np
import pytest
import h5py
from scipy.signal import StateSpace, lsim

from wecsim.environment.irregularWave import (
    jonswap_equal_energy_components, synthesize_irregular_response,
)
from wecsim.reference_models.wavestar import WaveStarLinkage


REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
HYDRO = os.environ.get("WEC_SIM_WAVESTAR_H5")
pytestmark = pytest.mark.skipif(
    not REFERENCE, reason="pinned WaveStar MATLAB output is absent",
)


def _source(name):
    path = Path(REFERENCE) / f"WECCCOMP_SOURCE_{name}.csv"
    data = np.loadtxt(path, delimiter=",", ndmin=2)
    assert np.isfinite(data).all()
    return data


def test_published_linkage_and_unforced_pto():
    float_body = _source("WECCCOMP_body1")
    arm = _source("WECCCOMP_body2")
    frame = _source("WECCCOMP_body3")
    pto = _source("WECCCOMP_pto1")
    assert all(item.shape == (14_121, 25)
               for item in (float_body, arm, frame, pto))
    np.testing.assert_allclose(float_body[:, 0], np.arange(14_121) * .01,
                               rtol=0, atol=1e-8)
    np.testing.assert_allclose(arm[:, 5], float_body[:, 5],
                               rtol=0, atol=1e-12)
    np.testing.assert_allclose(arm[:, 11], float_body[:, 11],
                               rtol=0, atol=1e-12)
    assert np.max(np.abs(frame[:, 1:13] - frame[0, 1:13])) < 1e-12

    linkage = WaveStarLinkage()
    angle = float_body[:, 5]
    angular_speed = float_body[:, 11]
    for response in (float_body, arm):
        neutral_center = response[0, [1, 3]]
        np.testing.assert_allclose(
            linkage.point_position(neutral_center, angle),
            response[:, [1, 3]], rtol=0, atol=1e-12,
        )
        np.testing.assert_allclose(
            linkage.point_velocity(neutral_center, angle, angular_speed),
            response[:, [7, 9]], rtol=0, atol=1e-12,
        )
    np.testing.assert_allclose(linkage.pto_stroke(angle), pto[:, 3],
                               rtol=0, atol=1e-12)
    np.testing.assert_allclose(linkage.pto_speed(angle, angular_speed),
                               pto[:, 9], rtol=0, atol=1e-12)
    np.testing.assert_allclose(pto[:, 13:25], 0, rtol=0, atol=1e-12)


@pytest.mark.skipif(not HYDRO, reason="pinned WaveStar HDF5 is absent")
def test_published_irregular_wave_and_excitation():
    source_components = _source("components")
    source_wave = _source("wave")
    source_body = _source("WECCCOMP_body1")
    assert source_components.shape == (500, 4)
    assert source_wave.shape == (14_121, 2)
    components = jonswap_equal_energy_components(
        HYDRO, significant_height=.0625, peak_period=1.412,
        directions=np.array([0.]), spreading=np.array([1.]),
        gamma=1, phase=source_components[:, 3, None],
    )
    np.testing.assert_allclose(components.omega, source_components[:, 0],
                               rtol=0, atol=1e-10)
    np.testing.assert_allclose(components.spectral_amplitude,
                               source_components[:, 1], rtol=0, atol=1e-12)
    np.testing.assert_allclose(components.d_omega, source_components[:, 2],
                               rtol=0, atol=1e-10)
    incident = synthesize_irregular_response(
        HYDRO, components, dt=.01, end_time=141.2,
        ramp_time=7.06, rho=1000, g=9.81,
    )
    np.testing.assert_allclose(incident.elevation, source_wave[:, 1],
                               rtol=0, atol=1e-10)
    np.testing.assert_allclose(incident.excitation_force,
                               source_body[:, 19:25], rtol=0, atol=1e-7)


@pytest.mark.skipif(not HYDRO, reason="pinned WaveStar HDF5 is absent")
def test_published_state_space_radiation_and_force_balance():
    body = _source("WECCCOMP_body1")
    forces = _source("body1_forces")
    assert forces.shape == (14_121, 37)
    with h5py.File(HYDRO) as h5:
        center = np.asarray(h5["body1/properties/cg"]).ravel()
        np.testing.assert_allclose(center[[0, 2]], body[0, [1, 3]],
                                   rtol=0, atol=1e-10)
        rho = float(np.asarray(h5["simulation_parameters/rho"]).item())
        fit = h5["body1/hydro_coeffs/radiation_damping/state_space"]
        order = np.asarray(fit["it"], dtype=int)
        matrix_a = np.asarray(fit["A/all"])
        matrix_b = np.asarray(fit["B/all"])
        matrix_c = np.asarray(fit["C/all"]) * rho
        frequency = np.r_[0, np.asarray(h5["simulation_parameters/w"]).ravel()]
        radiation = np.zeros((len(body), 6))
        for output in range(6):
            for input_axis in (0, 2, 4):
                count = order[output, input_axis]
                if count == 0:
                    continue
                state = StateSpace(
                    matrix_a[output, input_axis, :count, :count],
                    matrix_b[output, input_axis, :count, :1],
                    matrix_c[output, input_axis, :1, :count],
                    np.zeros((1, 1)),
                )
                _, signal, _ = lsim(
                    state, body[:, 7 + input_axis], body[:, 0], interp=True,
                )
                radiation[:, output] += signal
    joint_jacobian = np.array([
        body[0, 3] - .302, 0, -.438 - body[0, 1], 0, 1, 0,
    ])
    active_damping = np.zeros(len(frequency), dtype=complex)
    for output in (0, 2, 4):
        for input_axis in (0, 2, 4):
            count = order[output, input_axis]
            if count:
                state_a = matrix_a[output, input_axis, :count, :count]
                state_b = matrix_b[output, input_axis, :count, 0]
                state_c = matrix_c[output, input_axis, 0, :count]
                transfer = np.array([
                    state_c @ np.linalg.solve(
                        1j * omega * np.eye(count) - state_a, state_b,
                    ) for omega in frequency
                ])
                active_damping += (joint_jacobian[output]
                                   * joint_jacobian[input_axis] * transfer)
    assert active_damping.real.min() > 0
    np.testing.assert_allclose(radiation, forces[:, 1:7],
                               rtol=0, atol=.003)
    predicted_total = (body[:, 19:25] - forces[:, 1:7]
                       - forces[:, 7:13] - forces[:, 13:19]
                       - forces[:, 19:25] - forces[:, 25:31])
    np.testing.assert_allclose(predicted_total, body[:, 13:19],
                               rtol=0, atol=1e-9)

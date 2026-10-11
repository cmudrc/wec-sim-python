"""Isolate the published WaveStar ode8 step before its PTO activates."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim.irregularWave import jonswap_equal_energy_components
from wecsim.wavestar import run_wavestar_published


SOURCE = os.environ.get("WEC_SIM_WAVESTAR_ODE8_STEP_DIR")
HYDRO = os.environ.get("WEC_SIM_WAVESTAR_H5")
pytestmark = pytest.mark.skipif(
    not (SOURCE and HYDRO), reason="fixed-fit WaveStar ode8 step runs are absent",
)


def _read(folder, model, label):
    values = np.loadtxt(
        Path(SOURCE) / folder / f"{model}_{label}.csv",
        delimiter=",", ndmin=2,
    )
    assert np.isfinite(values).all()
    return values


def test_precontrol_ode8_step_audit():
    runs = []
    for folder, name, step in (
        ("50ms", "WECCCOMP_NMPC_ODE8_STEP", .05),
        ("25ms", "WECCCOMP_NMPC_ODE8_STEP", .025),
        ("12p5ms", "WECCCOMP_NMPC_ODE8_STEP", .0125),
    ):
        body = _read(folder, name, "WECCCOMP_Nonlinear_Model_Predictive_body1")
        pto = _read(folder, name, "WECCCOMP_Nonlinear_Model_Predictive_pto1")
        components = _read(folder, name, "components")
        wave = _read(folder, name, "wave")
        expected_count = round(9.95 / step) + 1
        assert body.shape == (expected_count, 25)
        assert pto.shape == (expected_count, 43)
        assert components.shape == (500, 4)
        assert wave.shape == (expected_count, 2)
        expected_time = np.arange(expected_count) * step
        for recorded in (body[:, 0], pto[:, 0], wave[:, 0]):
            np.testing.assert_allclose(recorded, expected_time,
                                       rtol=0, atol=1e-9)
        # Every source motion sample precedes the published 10 s start.
        np.testing.assert_allclose(pto[:, 31:37], 0., rtol=0, atol=1e-12)
        runs.append((body, components, wave, round(.05 / step)))

    reference_body, reference_components, reference_wave, _ = runs[0]
    for body, components, wave, stride in runs[1:]:
        np.testing.assert_allclose(components, reference_components,
                                   rtol=0, atol=1e-12)
        np.testing.assert_allclose(wave[::stride], reference_wave,
                                   rtol=0, atol=1e-10)
        np.testing.assert_allclose(body[::stride, 0], reference_body[:, 0],
                                   rtol=0, atol=1e-9)

    sea = jonswap_equal_energy_components(
        HYDRO, significant_height=.1042, peak_period=1.836,
        directions=np.array([0.]), spreading=np.array([1.]),
        gamma=3.3, phase=reference_components[:, 3, None],
    )
    python = run_wavestar_published(
        HYDRO, sea, dt=.001, end_time=9.95, ramp_time=25,
        g=9.80665, output_stride=50,
    )
    np.testing.assert_allclose(python.time, reference_body[:, 0],
                               rtol=0, atol=1e-9)
    np.testing.assert_allclose(python.wave_elevation, reference_wave[:, 1],
                               rtol=0, atol=1e-10)

    pitch = [body[::stride, 5] for body, _, _, stride in runs]
    speed = [body[::stride, 11] for body, _, _, stride in runs]
    errors = {
        "MATLAB_50_vs_25_ms_pitch_rad": float(np.max(np.abs(pitch[0] - pitch[1]))),
        "MATLAB_25_vs_12p5_ms_pitch_rad": float(np.max(np.abs(pitch[1] - pitch[2]))),
        "MATLAB_50_vs_25_ms_speed_rad_s": float(np.max(np.abs(speed[0] - speed[1]))),
        "MATLAB_25_vs_12p5_ms_speed_rad_s": float(np.max(np.abs(speed[1] - speed[2]))),
    }
    for folder, angle, angular_speed in zip(
        ("50ms", "25ms", "12p5ms"), pitch, speed,
    ):
        errors[f"Python_vs_MATLAB_{folder}_pitch_rad"] = float(
            np.max(np.abs(python.angle - angle))
        )
        errors[f"Python_vs_MATLAB_{folder}_speed_rad_s"] = float(
            np.max(np.abs(python.angular_speed - angular_speed))
        )
    assert errors["MATLAB_25_vs_12p5_ms_pitch_rad"] < .001
    assert (errors["MATLAB_25_vs_12p5_ms_pitch_rad"]
            < .5 * errors["MATLAB_50_vs_25_ms_pitch_rad"])
    assert (errors["MATLAB_25_vs_12p5_ms_speed_rad_s"]
            < .5 * errors["MATLAB_50_vs_25_ms_speed_rad_s"])
    for quantity, limit in (("pitch_rad", .0003),
                            ("speed_rad_s", .003)):
        base = errors[f"Python_vs_MATLAB_50ms_{quantity}"]
        half = errors[f"Python_vs_MATLAB_25ms_{quantity}"]
        quarter = errors[f"Python_vs_MATLAB_12p5ms_{quantity}"]
        assert quarter < limit
        assert quarter < .5 * half < .25 * base
    print("fixed-fit ode8 pre-control step audit:", errors)

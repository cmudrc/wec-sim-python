"""Isolate WaveStar's ode8 integration step before PTO control begins."""

import json
import os
from pathlib import Path

import numpy as np

from wecsim.dynamics.wavestar import run_wavestar_published
from wecsim.environment.irregularWave import jonswap_equal_energy_components


FIXTURE = Path(__file__).parent / "fixtures/matlab/wavestar_fixed_input_step"
SOURCE = os.environ.get("WEC_SIM_WAVESTAR_INTERNAL_STEP_DIR")
HYDRO = os.environ.get("WEC_SIM_WAVESTAR_H5", FIXTURE / "wavestar.h5")
PREFIX = "WECCCOMP_NMPC_INTERNAL_STEP_"


def _read(folder, label):
    if SOURCE:
        values = np.loadtxt(
            Path(SOURCE) / folder / f"{PREFIX}{label}.csv",
            delimiter=",", ndmin=2,
        )
    else:
        with np.load(FIXTURE / "source.npz") as saved:
            values = saved[f"{folder}_{label}"]
    assert np.isfinite(values).all(), (folder, label)
    return values


def test_fixed_inputs_isolate_precontrol_ode8_substeps():
    runs = []
    common_time = np.arange(200) * .05
    for folder, solver_step in (("50ms", .05), ("25ms", .025),
                                ("12p5ms", .0125)):
        body = _read(folder, "WECCCOMP_Nonlinear_Model_Predictive_body1")
        pto = _read(folder, "WECCCOMP_Nonlinear_Model_Predictive_pto1")
        components = _read(folder, "components")
        wave = _read(folder, "wave")
        steps = _read(folder, "steps")
        assert body.shape == (200, 25)
        # The PTO block logs at Simulink's internal step, while WEC-Sim body
        # and wave outputs retain the published 50 ms reporting interval.
        pto_time = np.arange(round(9.95 / solver_step) + 1) * solver_step
        assert pto.shape == (pto_time.size, 43)
        assert components.shape == (500, 4)
        assert wave.shape == (200, 2)
        np.testing.assert_allclose(steps, [[.05, .05, solver_step, solver_step]],
                                   rtol=0, atol=1e-12)
        for recorded in (body[:, 0], wave[:, 0]):
            np.testing.assert_allclose(recorded, common_time,
                                       rtol=0, atol=1e-9)
        np.testing.assert_allclose(pto[:, 0], pto_time, rtol=0, atol=1e-9)
        np.testing.assert_allclose(pto[:, 31:37], 0., rtol=0, atol=1e-12)
        runs.append((body, components, wave))

    reference_body, reference_components, reference_wave = runs[0]
    for _, components, wave in runs[1:]:
        np.testing.assert_allclose(components, reference_components,
                                   rtol=0, atol=1e-12)
        np.testing.assert_allclose(wave, reference_wave,
                                   rtol=0, atol=1e-12)

    sea = jonswap_equal_energy_components(
        HYDRO, significant_height=.1042, peak_period=1.836,
        directions=np.array([0.]), spreading=np.array([1.]),
        gamma=3.3, phase=reference_components[:, 3, None],
    )
    python = run_wavestar_published(
        HYDRO, sea, dt=.001, end_time=9.95, ramp_time=25,
        g=9.80665, output_stride=50,
    )
    np.testing.assert_allclose(python.time, common_time, rtol=0, atol=1e-9)
    np.testing.assert_allclose(python.wave_elevation, reference_wave[:, 1],
                               rtol=0, atol=1e-10)

    pitch = [body[:, 5] for body, _, _ in runs]
    speed = [body[:, 11] for body, _, _ in runs]
    errors = {
        "source_50_to_25_pitch_rad": float(np.max(np.abs(pitch[0] - pitch[1]))),
        "source_25_to_12p5_pitch_rad": float(np.max(np.abs(pitch[1] - pitch[2]))),
        "source_50_to_25_speed_rad_s": float(np.max(np.abs(speed[0] - speed[1]))),
        "source_25_to_12p5_speed_rad_s": float(np.max(np.abs(speed[1] - speed[2]))),
    }
    for folder, angle, angular_speed in zip(
        ("50ms", "25ms", "12p5ms"), pitch, speed,
    ):
        errors[f"python_vs_{folder}_pitch_rad"] = float(
            np.max(np.abs(python.angle - angle))
        )
        errors[f"python_vs_{folder}_speed_rad_s"] = float(
            np.max(np.abs(python.angular_speed - angular_speed))
        )
    # The published-step gap is an identity guard for this archived input.
    # The second two limits were declared before the hosted source run: they
    # test whether internal integration alone explains most of that gap.
    assert .0025 < errors["python_vs_50ms_pitch_rad"] < .0045
    assert np.max(np.abs(pitch[0] - pitch[2])) > .0017
    assert errors["python_vs_12p5ms_pitch_rad"] < .001
    print("fixed-input WaveStar solver-step audit:", json.dumps(errors,
                                                             sort_keys=True))

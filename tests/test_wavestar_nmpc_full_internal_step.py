"""Separate WaveStar's internal ode8 step from its 50 ms NMPC schedule."""

import json
import os
from pathlib import Path

import h5py
import numpy as np
import pytest
from scipy.io import loadmat

from wecsim import WaveStarNmpcPTO
from wecsim.dynamics.wavestar import run_wavestar_published
from wecsim.environment.irregularWave import jonswap_equal_energy_components


FIXTURE = Path(__file__).parent / "fixtures/matlab/wavestar_nmpc_full_step"
PRECONTROL = Path(__file__).parent / "fixtures/matlab/wavestar_fixed_input_step"
SOURCE = os.environ.get("WEC_SIM_WAVESTAR_FULL_INTERNAL_STEP_DIR")
HYDRO = os.environ.get("WEC_SIM_WAVESTAR_H5")
PREFIX = "WECCCOMP_NMPC_INTERNAL_STEP_FULL_"


def _fresh(label):
    values = np.loadtxt(Path(SOURCE) / f"{PREFIX}{label}.csv",
                        delimiter=",", ndmin=2)
    assert np.isfinite(values).all(), label
    return values


def _difference_metrics(python, source, count):
    difference = np.asarray(python[:count]) - np.asarray(source[:count])
    return {"maximum": float(np.max(np.abs(difference))),
            "rms": float(np.sqrt(np.mean(difference**2)))}


def test_archived_full_baseline_reuses_the_fixed_precontrol_input():
    with np.load(FIXTURE / "baseline.npz") as archived, \
            np.load(PRECONTROL / "source.npz") as short:
        assert archived["body"].shape == (4501, 25)
        assert archived["pto"].shape == (4501, 43)
        assert archived["components"].shape == (500, 4)
        np.testing.assert_array_equal(
            archived["body"][:200],
            short["50ms_WECCCOMP_Nonlinear_Model_Predictive_body1"],
        )
        np.testing.assert_array_equal(
            archived["pto"][:200],
            short["50ms_WECCCOMP_Nonlinear_Model_Predictive_pto1"],
        )
        np.testing.assert_array_equal(
            archived["components"], short["50ms_components"],
        )
        np.testing.assert_array_equal(
            archived["wave"][:200], short["50ms_wave"],
        )
        for label in ("command", "power", "energy"):
            np.testing.assert_allclose(
                archived[f"{label}_time"], np.arange(4501) * .05,
                rtol=0, atol=1e-9,
            )


@pytest.mark.skipif(not (SOURCE and HYDRO),
                    reason="full-duration 25 ms MATLAB source run absent")
def test_full_internal_step_with_published_controller_grid():
    with np.load(FIXTURE / "baseline.npz") as archived:
        baseline = {key: archived[key] for key in archived.files}
    with np.load(PRECONTROL / "source.npz") as saved:
        short = {key: saved[key] for key in saved.files}
    body = _fresh("WECCCOMP_Nonlinear_Model_Predictive_body1")
    pto = _fresh("WECCCOMP_Nonlinear_Model_Predictive_pto1")
    components = _fresh("components")
    wave = _fresh("wave")
    steps = _fresh("steps")
    controller = loadmat(Path(SOURCE) / f"{PREFIX}controller.mat",
                         simplify_cells=True)
    command = np.asarray(controller["cmd_ptoM"]["signals"]["values"])
    command_time = np.asarray(controller["cmd_ptoM"]["time"])
    source_power = np.asarray(
        controller["Output_power"]["signals"]["values"]
    )
    source_energy = np.asarray(
        controller["Output_energy"]["signals"]["values"]
    )
    assert body.shape == (4501, 25)
    assert pto.shape == (9001, 43)
    assert wave.shape == (4501, 2)
    assert command.shape == command_time.shape == (4501,)
    assert source_power.shape == source_energy.shape == (4501,)
    np.testing.assert_allclose(steps, [[.05, .05, .025, .025]],
                               rtol=0, atol=1e-12)
    np.testing.assert_allclose(body[:, 0], np.arange(4501) * .05,
                               rtol=0, atol=1e-9)
    np.testing.assert_allclose(pto[:, 0], np.arange(9001) * .025,
                               rtol=0, atol=1e-9)
    np.testing.assert_allclose(command_time, body[:, 0], rtol=0, atol=1e-9)
    for signal in ("Output_power", "Output_energy"):
        np.testing.assert_allclose(controller[signal]["time"], command_time,
                                   rtol=0, atol=1e-9)
    np.testing.assert_allclose(wave, baseline["wave"], rtol=0, atol=1e-12)
    np.testing.assert_allclose(components, baseline["components"],
                               rtol=0, atol=1e-12)
    np.testing.assert_allclose(
        body[:200], short["25ms_WECCCOMP_Nonlinear_Model_Predictive_body1"],
        rtol=0, atol=1e-9,
    )
    np.testing.assert_allclose(
        pto[:399], short["25ms_WECCCOMP_Nonlinear_Model_Predictive_pto1"],
        rtol=0, atol=1e-9,
    )
    assert np.max(np.abs(body[:200, 5] - baseline["body"][:200, 5])) > .0017
    assert np.max(np.abs(command[command_time >= 15])) > 1

    # The baseline and fresh BEMIO containers may have different metadata;
    # compare every numerical input dataset rather than their file hashes.
    with h5py.File(HYDRO) as fresh, \
            h5py.File(PRECONTROL / "wavestar.h5") as pinned:
        fresh_keys, pinned_keys = [], []
        fresh.visititems(lambda key, value: fresh_keys.append(key)
                         if isinstance(value, h5py.Dataset) else None)
        pinned.visititems(lambda key, value: pinned_keys.append(key)
                          if isinstance(value, h5py.Dataset) else None)
        assert fresh_keys == pinned_keys
        for key in pinned_keys:
            np.testing.assert_array_equal(fresh[key][()], pinned[key][()])

    sea = jonswap_equal_energy_components(
        HYDRO, significant_height=.1042, peak_period=1.836,
        directions=np.array([0.]), spreading=np.array([1.]),
        gamma=3.3, phase=components[:, 3, None],
    )
    python_pto = WaveStarNmpcPTO(plant_dt=.001, control_dt=.05)
    python = run_wavestar_published(
        HYDRO, sea, dt=.001, end_time=225, ramp_time=25,
        g=9.80665, pto_controller=python_pto, output_stride=50,
    )
    np.testing.assert_allclose(python.time, body[:, 0], rtol=0, atol=1e-9)
    np.testing.assert_allclose(python_pto.command_time, command_time,
                               rtol=0, atol=1e-9)
    assert np.isfinite(python.angle).all()
    assert np.isfinite(python_pto.command_torque).all()

    sampled_pto = pto[::2]
    metrics = {}
    for label, count in (("precontrol", 200), ("first_30s", 601),
                         ("full_225s", 4501)):
        for channel, old, new in (
            ("pitch_rad", baseline["body"][:, 5], body[:, 5]),
            ("speed_rad_s", baseline["body"][:, 11], body[:, 11]),
            ("pto_force_N", baseline["pto"][:, 33], sampled_pto[:, 33]),
            ("command_Nm", baseline["command"], command),
            ("reported_power", baseline["power"], source_power),
        ):
            metrics[f"{label}_source_self_change_{channel}"] = (
                _difference_metrics(old, new, count)
            )
        for source_name, source_body, source_pto, source_command in (
            ("published_50ms", baseline["body"], baseline["pto"],
             baseline["command"]),
            ("internal_25ms", body, sampled_pto, command),
        ):
            for channel, actual, expected in (
                ("pitch_rad", python.angle, source_body[:, 5]),
                ("speed_rad_s", python.angular_speed, source_body[:, 11]),
                ("pto_force_N", python.pto_force, source_pto[:, 33]),
                ("command_Nm", python_pto.command_torque, source_command),
            ):
                metrics[f"{label}_{source_name}_{channel}"] = (
                    _difference_metrics(actual, expected, count)
                )
    for name, source_pto in (("published_50ms", baseline["pto"]),
                             ("internal_25ms", sampled_pto)):
        source_work = np.trapezoid(
            source_pto[:, 33] * source_pto[:, 9], x=source_pto[:, 0],
        )
        python_work = np.trapezoid(
            python.pto_force * python.pto_speed, x=python.time,
        )
        metrics[f"{name}_relative_mechanical_work"] = float(
            abs(python_work - source_work) / max(abs(source_work), 1e-12)
        )
    metrics["source_reported_energy_self_change_relative"] = float(
        abs(source_energy[-1] - baseline["energy"][-1])
        / max(abs(baseline["energy"][-1]), 1e-12)
    )
    assert .02 < metrics["full_225s_published_50ms_pitch_rad"]["maximum"] < .06
    print("WaveStar fixed-input 225 s internal-step audit:",
          json.dumps(metrics, sort_keys=True))

"""Diagnose WaveStar plant-step effects with its controller clock held fixed."""

import json
import os
from pathlib import Path

import numpy as np
import pytest

from wecsim import WaveStarNmpcPTO
from wecsim.dynamics.wavestar import run_wavestar_published
from wecsim.environment.irregularWave import jonswap_equal_energy_components


ROOT = Path(__file__).parent
SOURCE = os.environ.get("WEC_SIM_WAVESTAR_HELD_CLOCK_DIR")
BASELINE = ROOT / "fixtures/matlab/wavestar_nmpc_full_step/baseline.npz"
HYDRO = ROOT / "fixtures/matlab/wavestar_fixed_input_step/wavestar.h5"


def _source(name):
    values = np.loadtxt(Path(SOURCE) / f"{name}.csv", delimiter=",", ndmin=2)
    assert np.isfinite(values).all(), name
    return values


def _gap(first, second):
    difference = np.asarray(first) - np.asarray(second)
    return {"maximum": float(np.max(np.abs(difference))),
            "rms": float(np.sqrt(np.mean(difference**2)))}


@pytest.mark.skipif(not SOURCE, reason="fresh held-clock MATLAB run absent")
def test_full_held_clock_measures_source_and_python_differences():
    with np.load(BASELINE, allow_pickle=False) as saved:
        published = {name: saved[name].copy() for name in saved.files}
    body = _source("float")
    axial = _source("pto_axial")
    command = _source("command")
    wave = _source("wave")
    components = _source("components")
    assert body.shape == (4501, 13)
    assert axial.shape in ((4501, 3), (9001, 3))
    assert command.shape in ((4501, 2), (9001, 2))
    assert wave.shape == (4501, 2)
    assert components.shape == (500, 4)
    np.testing.assert_allclose(body[:, 0], np.arange(4501) * .05,
                               rtol=0, atol=1e-9)
    np.testing.assert_allclose(axial[:, 0], np.linspace(0, 225, len(axial)),
                               rtol=0, atol=1e-9)
    np.testing.assert_allclose(command[:, 0], np.linspace(0, 225, len(command)),
                               rtol=0, atol=1e-9)
    np.testing.assert_allclose(wave, published["wave"], rtol=0, atol=1e-12)
    np.testing.assert_allclose(components, published["components"],
                               rtol=0, atol=1e-12)
    if len(command) == 9001:
        np.testing.assert_allclose(command[1::2, 1], command[:-1:2, 1],
                                   rtol=0, atol=1e-9)
        held_command = command[::2, 1]
    else:
        held_command = command[:, 1]
    assert np.ptp(held_command[601:]) > 1  # The controller remains active.
    assert np.max(np.abs(axial[:, 1])) > 1

    sea = jonswap_equal_energy_components(
        HYDRO, significant_height=.1042, peak_period=1.836,
        directions=np.array([0.]), spreading=np.array([1.]), gamma=3.3,
        phase=components[:, 3, None],
    )
    python = run_wavestar_published(
        HYDRO, sea, dt=.001, end_time=225, ramp_time=25, g=9.80665,
        pto_controller=WaveStarNmpcPTO(plant_dt=.001, control_dt=.05),
        output_stride=50,
    )
    np.testing.assert_allclose(python.time, body[:, 0], rtol=0, atol=1e-9)
    axial_on_output = axial[::2] if len(axial) == 9001 else axial
    metrics = {}
    for window, count in (("first_10s", 200), ("first_30s", 601),
                          ("full_225s", 4501)):
        for name, old, held, independent in (
            ("pitch_rad", published["body"][:, 5], body[:, 5], python.angle),
            ("pitch_speed_rad_s", published["body"][:, 11], body[:, 11],
             python.angular_speed),
            ("pto_force_N", published["pto"][:, 33], axial_on_output[:, 1],
             python.pto_force),
        ):
            metrics[f"{window}_{name}_source_self_change"] = _gap(
                old[:count], held[:count])
            metrics[f"{window}_{name}_python_vs_held"] = _gap(
                independent[:count], held[:count])
        metrics[f"{window}_command_source_self_change"] = _gap(
            published["command"][:count], held_command[:count])
    old_work = np.trapezoid(published["pto"][:, 33]
                            * published["pto"][:, 9],
                            x=published["pto"][:, 0])
    held_work = np.trapezoid(axial[:, 1] * axial[:, 2], x=axial[:, 0])
    python_work = np.trapezoid(python.pto_force * python.pto_speed,
                               x=python.time)
    metrics["relative_pto_work_source_self_change"] = float(
        abs(held_work - old_work) / max(abs(old_work), 1e-12))
    metrics["relative_pto_work_python_vs_held"] = float(
        abs(python_work - held_work) / max(abs(held_work), 1e-12))
    assert all(np.isfinite(value) for entry in metrics.values()
               for value in (entry.values() if isinstance(entry, dict)
                             else (entry,)))
    print("WaveStar full held-clock diagnostic:", json.dumps(metrics, sort_keys=True))

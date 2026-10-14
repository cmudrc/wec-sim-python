"""Independent WaveBot motion against the pinned 600 s MATLAB application."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim import run_wavebot_impedance


REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
HYDRO = os.environ.get("WEC_SIM_WAVEBOT_H5")
MULTISINE = os.environ.get("WEC_SIM_WAVEBOT_MULTISINE")
PORTABLE = Path(__file__).parent / "fixtures/matlab/wavebot_impedance"


def _source(name):
    values = np.loadtxt(
        Path(REFERENCE) / f"WAVEBOT_IMPEDANCE_{name}.csv",
        delimiter=",", ndmin=2,
    )
    assert np.isfinite(values).all(), name
    return values


@pytest.fixture(scope="module")
def comparison():
    if any((REFERENCE, HYDRO, MULTISINE)) and not all((REFERENCE, HYDRO, MULTISINE)):
        pytest.fail("live WaveBot comparison needs source output, HDF5, and multisine")
    if REFERENCE:
        body = _source("CalcImpedance_body1")
        force = _source("body1_forces")
        mooring = _source("mooring")
        actuation = _source("actuation")
        assert body.shape == (60_001, 25)
        assert force.shape == (60_001, 37)
        assert mooring.shape == (60_001, 19)
        reference = {
            "body_time": body[:, 0],
            "body_position": body[:, [1, 3, 5]],
            "body_velocity": body[:, [7, 9, 11]],
            "actuation": actuation,
            "mooring_force": mooring[:, [13, 15, 17]],
            "linear_damping_force": -force[:, [25, 27, 29]],
            "drag_force": -force[:, [19, 21, 23]],
        }
        hydro, multisine = HYDRO, MULTISINE
    else:
        with np.load(PORTABLE / "published_600s.npz") as saved:
            reference = {name: saved[name].copy() for name in saved.files}
        hydro = PORTABLE / "waveBotBuoy.h5"
        multisine = PORTABLE / "multisine3DOFA.mat"
    return reference, run_wavebot_impedance(
        hydro, multisine, source_linear_damping=True,
    ), hydro, multisine


def test_published_wavebot_coupled_motion_and_force_channels(comparison):
    source, result, _, _ = comparison
    assert source["body_time"].shape == (60_001,)
    np.testing.assert_allclose(result.time, source["body_time"], rtol=0, atol=1e-8)

    active = [0, 2, 4]
    position_error = np.max(
        np.abs(result.body_position[:, active] - source["body_position"]), axis=0,
    )
    velocity_error = np.max(
        np.abs(result.body_velocity[:, active] - source["body_velocity"]), axis=0,
    )
    np.testing.assert_array_less(position_error, [.0002, .0003, .0004])
    np.testing.assert_array_less(velocity_error, [.0012, .001, .0018])

    sampled_command = np.column_stack([
        np.interp(result.time, source["actuation"][:, 0],
                  source["actuation"][:, i])
        for i in (1, 2, 3)
    ])
    sampled_command[:, 2] *= -1
    np.testing.assert_allclose(result.applied_force, sampled_command,
                               rtol=0, atol=1e-8)
    force_pairs = (
        (result.mooring_force, source["mooring_force"], [4, 1.3, .4]),
        (result.linear_damping_force,
         source["linear_damping_force"], [1.3, 1.3, .15]),
        (result.drag_force, source["drag_force"], [.2, .8, .5]),
    )
    for python_force, matlab_force, limits in force_pairs:
        error = np.max(np.abs(python_force[:, active] - matlab_force), axis=0)
        np.testing.assert_array_less(error, limits)


def test_ordinary_diagonal_damping_stays_the_default(comparison):
    source, source_compatible, hydro, multisine = comparison
    physical = run_wavebot_impedance(hydro, multisine)
    # Each positive diagonal coefficient dissipates power independently.
    damping_power = np.sum(
        physical.linear_damping_force * physical.body_velocity, axis=1,
    )
    assert np.max(damping_power) < 1e-10
    # MATLAB's linearly indexed matrix is different from that physical
    # default. Keep the source convention isolated and report both traces.
    assert np.max(np.abs(physical.body_position[:, 2]
                         - source["body_position"][:, 1])) > .01
    assert np.max(np.abs(physical.body_position[:, 4]
                         - source["body_position"][:, 2])) > .005
    assert np.max(np.abs(source_compatible.body_position[:, 2]
                         - source["body_position"][:, 1])) < .0003

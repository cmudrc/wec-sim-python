"""Keep the published three-body MBARI cable trajectory paired in Python CI."""

from pathlib import Path

import numpy as np

from wecsim import run_mbari_cable


ROOT = Path(__file__).resolve().parents[1]
HYDRO = ROOT / "examples/data/mbari.h5"
MATLAB = ROOT / "tests/fixtures/matlab/mbari_cable.npz"


def test_published_mbari_coupled_motion_and_cable_force():
    """Advance from the pinned HDF5, without source motion or force inputs."""
    result = run_mbari_cable(HYDRO)
    with np.load(MATLAB, allow_pickle=False) as source:
        assert source["time"].shape == (12_001,)
        assert source["body_position"].shape == (12_001, 3, 6)
        assert source["body_velocity"].shape == (12_001, 3, 6)
        np.testing.assert_allclose(result.time, source["time"],
                                   rtol=0, atol=1e-10)
        np.testing.assert_allclose(result.wave_elevation,
                                   source["wave_elevation"],
                                   rtol=0, atol=1e-11)

        position_error = np.max(np.abs(
            result.body_position - source["body_position"],
        ), axis=0)
        velocity_error = np.max(np.abs(
            result.body_velocity - source["body_velocity"],
        ), axis=0)
        # These are the original full-duration gates, including inactive
        # axes and the brief, stiff cable-snap velocity transient.
        position_gate = np.array([
            [.01, 1e-10, .05, 1e-10, .01, 1e-10],
            [.005, 1e-10, .05, 1e-10, .003, 1e-10],
            [.01, 1e-10, .05, 1e-10, .003, 1e-10],
        ])
        velocity_gate = np.array([
            [.03, 1e-10, .45, 1e-10, .05, 1e-10],
            [.015, 1e-10, .35, 1e-10, .01, 1e-10],
            [.04, 1e-10, .45, 1e-10, .02, 1e-10],
        ])
        np.testing.assert_array_less(position_error, position_gate)
        np.testing.assert_array_less(velocity_error, velocity_gate)

        displacement_error = (result.cable_displacement
                              - source["cable_displacement"])
        speed_error = result.cable_speed - source["cable_speed"]
        force_error = result.cable_force_z - source["cable_force_z"]
        assert np.max(np.abs(displacement_error)) < .06
        assert np.max(np.abs(speed_error)) < .8
        assert np.max(np.abs(force_error)) < 45_000
        assert np.sqrt(np.mean(force_error**2)) < 3_500
        assert abs(np.trapezoid(force_error, result.time)) < 500
        assert abs(np.mean(result.cable_force_z != 0)
                   - np.mean(source["cable_force_z"] != 0)) < .02

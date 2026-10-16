"""Isolate the published RM3 MoorDyn solver's startup-step sensitivity."""

import os
from pathlib import Path
import shutil

import numpy as np
import pytest

from wecsim import JONSWAPWave, MoorDyn, WEC


FIXTURE = (Path(__file__).parent / "fixtures/matlab"
           / "rm3_moordyn_initial_step_2s/signals.npz")
PUBLISHED = (Path(__file__).parent / "fixtures/matlab"
             / "rm3_moordyn_source_step.npz")


def test_initial_step_only_changes_source_surge_and_pitch():
    with np.load(FIXTURE, allow_pickle=False) as source:
        time = source["time"]
        default = source["default_position"]
        fixed = source["initial_step_position"]
        assert default.shape == fixed.shape == (2, 201, 6)
        assert source["default_velocity"].shape == (2, 201, 6)
        assert source["initial_step_velocity"].shape == (2, 201, 6)
        np.testing.assert_allclose(time, np.arange(201) * .01,
                                   rtol=0, atol=1e-14)
        np.testing.assert_array_equal(source["default_phase"],
                                      source["initial_step_phase"])
        np.testing.assert_array_equal(source["default_wave"],
                                      source["initial_step_wave"])
        assert source["default_phase"].shape == (1000,)
        with np.load(PUBLISHED, allow_pickle=False) as published:
            np.testing.assert_array_equal(source["default_phase"],
                                          published["phase"][0])
            np.testing.assert_allclose(source["default_wave"][::10],
                                       published["wave"][0, :21],
                                       rtol=0, atol=2e-14)
            for body in (0, 1):
                np.testing.assert_allclose(
                    default[body, ::10],
                    published[f"body{body + 1}_position"][0, :21],
                    rtol=0, atol=1e-12,
                )

    # This locks down a controlled source numerical result, not a physical
    # accuracy or published-input parity claim.
    surge_shift = np.max(np.abs(fixed[:, :, 0] - default[:, :, 0]), axis=1)
    pitch_shift = np.max(np.abs(fixed[:, :, 4] - default[:, :, 4]))
    np.testing.assert_allclose(surge_shift, [.036418798365, .062001492983],
                               rtol=0, atol=1e-9)
    assert .0012 < pitch_shift < .0013


@pytest.mark.skipif(
    not all(os.getenv(name) for name in (
        "WEC_SIM_RM3_H5", "WEC_SIM_MOORDYN_LIBRARY",
        "WEC_SIM_RM3_MOORDYN_LINES",
    )),
    reason="pinned RM3 HDF5, native MoorDyn, and lines input not provided",
)
def test_independent_python_motion_after_source_initial_step_change(tmp_path):
    """Run live mooring feedback against both source solver settings."""
    with np.load(FIXTURE, allow_pickle=False) as source:
        phase = source["default_phase"]
        wave = source["default_wave"]
        original = source["default_position"]
        changed = source["initial_step_position"]
    phase_path = tmp_path / "phase.csv"
    np.savetxt(phase_path, phase, delimiter=",")
    line_dir = tmp_path / "Mooring"
    line_dir.mkdir()
    lines = line_dir / "lines.txt"
    shutil.copyfile(os.environ["WEC_SIM_RM3_MOORDYN_LINES"], lines)
    hydro = Path(os.environ["WEC_SIM_RM3_H5"])
    wec = WEC("RM3 initial-step diagnostic")
    float_body = wec.body("float", hydro, inertia=(0, 21_306_090.66, 0))
    spar = wec.body("spar", hydro, inertia=(0, 94_407_091.24, 0))
    wec.floating_joint(
        float_body, spar, damping=1_200_000,
        moordyn=MoorDyn(os.environ["WEC_SIM_MOORDYN_LIBRARY"], lines),
        moordyn_point=spar.at(0, 0, 21.5),
    )
    result = wec.run(
        JONSWAPWave(2, 8, phase_file=phase_path,
                    discretization="traditional"),
        dt=.01, end_time=2, ramp_time=0, radiation_memory=60,
        initial_coordinate={"spar_heave": -.21},
    )
    np.testing.assert_allclose(result.wave_elevation, wave,
                               rtol=0, atol=2e-13)
    for body, name in enumerate(("float", "spar")):
        pose = result.bodies[name].position
        surge_error = np.max(np.abs(pose[:, 0] - changed[body, :, 0]))
        original_error = np.max(np.abs(pose[:, 0] - original[body, :, 0]))
        pitch_error = np.max(np.abs(pose[:, 4] - changed[body, :, 4]))
        print(f"{name}: initial-step surge={surge_error:.8g} m, "
              f"default surge={original_error:.8g} m, "
              f"initial-step pitch={pitch_error:.8g} rad")
        assert surge_error < .003
        assert surge_error < .1 * original_error
        assert pitch_error < 5e-5

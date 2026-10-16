"""Keep the full RM3 MoorDyn ode45 startup-step diagnostic reproducible."""

import os
from pathlib import Path
import shutil

import numpy as np
import pytest

from wecsim import JONSWAPWave, MoorDyn, WEC


ROOT = Path(__file__).parent
FIXTURE = ROOT / "fixtures/matlab/rm3_moordyn_initial_step_80s/signals.npz"
PUBLISHED = ROOT / "fixtures/matlab/rm3_moordyn_source_step.npz"


def _maximum(left, right):
    assert left.shape == right.shape
    return float(np.max(np.abs(left - right)))


def _squared_gap_reduction(published, changed, refined):
    original = refined - published
    assert np.sum(original**2) > 0
    return 1 - np.sum((refined - changed)**2) / np.sum(original**2)


def test_initial_step_explains_full_duration_source_surge_and_pitch_gap():
    """The unchanged local arm must first reproduce the archived source."""
    with (np.load(FIXTURE, allow_pickle=False) as run,
          np.load(PUBLISHED, allow_pickle=False) as archive):
        np.testing.assert_allclose(run["time"], np.arange(801) * .1,
                                   rtol=0, atol=2e-14)
        np.testing.assert_allclose(run["fairlead_time"],
                                   np.arange(8001) * .01,
                                   rtol=0, atol=2e-12)
        np.testing.assert_array_equal(run["phase"][0], run["phase"][1])
        np.testing.assert_array_equal(run["phase"][0], archive["phase"][0])
        np.testing.assert_array_equal(run["wave"][0], run["wave"][1])
        np.testing.assert_allclose(run["wave"][0, :, 1],
                                   archive["wave"][0], rtol=0, atol=1e-13)
        assert run["body_position"].shape == (2, 2, 801, 6)
        assert run["body_velocity"].shape == (2, 2, 801, 6)
        assert run["mooring_force"].shape == (2, 801, 6)
        assert run["pto_force"].shape == (2, 801)
        assert run["fairlead_tension"].shape == (2, 8001, 3)

        for body in (0, 1):
            name = f"body{body + 1}_position"
            # This identifies the local auto-step arm as the original 80 s
            # trajectory; a short matching prefix would be weaker evidence.
            assert _maximum(run["body_position"][0, body],
                            archive[name][0]) < 3e-12
            published, changed, refined = (archive[name][0],
                                           run["body_position"][1, body],
                                           archive[name][2])
            assert _maximum(published[:, 0], changed[:, 0]) > .18
            assert _maximum(changed[:, 0], refined[:, 0]) < .008
            assert _squared_gap_reduction(published[:, 0], changed[:, 0],
                                          refined[:, 0]) > .999
            assert _maximum(published[:, 4], changed[:, 4]) > .0017
            assert _maximum(changed[:, 4], refined[:, 4]) < 5e-5
            assert _squared_gap_reduction(published[:, 4], changed[:, 4],
                                          refined[:, 4]) > .999

        assert _maximum(run["mooring_force"][0],
                        archive["mooring_force"][0]) < 1e-5
        published = archive["mooring_force"][0, :, 4]
        changed = run["mooring_force"][1, :, 4]
        refined = archive["mooring_force"][2, :, 4]
        assert _maximum(published, changed) > 100_000
        assert _maximum(changed, refined) < 10_000
        assert _squared_gap_reduction(published, changed, refined) > .995
        # InitialStep is not an all-channel explanation: spar heave has a
        # much smaller squared-gap reduction than surge and pitch.
        published = archive["body2_position"][0, :, 2]
        changed = run["body_position"][1, 1, :, 2]
        refined = archive["body2_position"][2, :, 2]
        assert _squared_gap_reduction(published, changed, refined) < .9


@pytest.mark.skipif(
    not all(os.getenv(name) for name in (
        "WEC_SIM_RM3_H5", "WEC_SIM_MOORDYN_LIBRARY",
        "WEC_SIM_RM3_MOORDYN_LINES",
    )),
    reason="pinned RM3 HDF5, native MoorDyn, and lines input not provided",
)
def test_independent_python_against_changed_80s_source(tmp_path):
    """Advance live mooring feedback; source motion supplies no solver input."""
    with np.load(FIXTURE, allow_pickle=False) as source:
        phase = source["phase"][0].copy()
        expected_wave = source["wave"][1, :, 1].copy()
        expected_position = source["body_position"][1].copy()
        expected_velocity = source["body_velocity"][1].copy()
        expected_mooring = source["mooring_force"][1].copy()
        expected_pto = source["pto_force"][1].copy()
        expected_fairlead = source["fairlead_tension"][1].copy()
    phase_path = tmp_path / "phase.csv"
    np.savetxt(phase_path, phase, delimiter=",")
    lines = tmp_path / "Mooring/lines.txt"
    lines.parent.mkdir()
    shutil.copyfile(os.environ["WEC_SIM_RM3_MOORDYN_LINES"], lines)
    wec = WEC("RM3 initial-step 80 s diagnostic")
    float_body = wec.body("float", os.environ["WEC_SIM_RM3_H5"],
                          inertia=(0, 21_306_090.66, 0))
    spar = wec.body("spar", os.environ["WEC_SIM_RM3_H5"],
                    inertia=(0, 94_407_091.24, 0))
    wec.floating_joint(
        float_body, spar, damping=1_200_000,
        moordyn=MoorDyn(os.environ["WEC_SIM_MOORDYN_LIBRARY"], lines),
        moordyn_point=spar.at(0, 0, 21.5),
    )
    result = wec.run(
        JONSWAPWave(2, 8, phase_file=phase_path,
                    discretization="traditional"),
        dt=.01, end_time=80, ramp_time=0, radiation_memory=60,
        initial_coordinate={"spar_heave": -.21},
    )
    indices = np.arange(801) * 10
    assert _maximum(result.wave_elevation[indices], expected_wave) < 2e-13
    for body, name in enumerate(("float", "spar")):
        for field, expected, limits in (
            ("position", expected_position[body], (.01, .0005, 7e-5)),
            ("velocity", expected_velocity[body], (.002, .0008, 5e-5)),
        ):
            actual = getattr(result.bodies[name], field)[indices]
            for axis, limit in zip((0, 2, 4), limits, strict=True):
                assert _maximum(actual[:, axis], expected[:, axis]) < limit, (
                    name, field, axis)
    force = dict(result.raw.extra_outputs)["moordyn_connection_force"][indices]
    for axis, limit in ((0, 1500), (2, 800), (4, 10_000)):
        assert _maximum(force[:, axis], expected_mooring[:, axis]) < limit
    assert _maximum(result.ptos["relative_heave"].force[indices],
                    expected_pto) < 800
    fairlead = np.loadtxt(lines.with_suffix(".out"), skiprows=1, ndmin=2)
    assert fairlead.shape == (8001, 4)
    for line, limit in enumerate((2500, 900, 900), start=1):
        assert _maximum(fairlead[:, line],
                        expected_fairlead[:, line - 1]) < limit

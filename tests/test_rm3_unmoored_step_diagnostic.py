"""Separate source solver sensitivity from MoorDyn in a derived RM3 sea."""

from pathlib import Path
from hashlib import sha256

import numpy as np

from wecsim import JONSWAPWave, WEC


ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests/fixtures/matlab/rm3_unmoored_step_10s"
SOURCE = FIXTURE / "source.npz"
HYDRO = FIXTURE / "rm3.h5"
HASHES = {
    SOURCE: "2ce203b8c6a6b90d1d2439618f71ab71b57eda3eba13e6b04afdff4998ea6d02",
    HYDRO: "a909360c6479773c113f792be72afb1351cae9299d62e501fa80561dca450aa3",
}


def test_unmoored_rm3_source_step_and_independent_python_motion(tmp_path):
    for path, digest in HASHES.items():
        assert sha256(path.read_bytes()).hexdigest() == digest
    with np.load(SOURCE, allow_pickle=False) as fixture:
        source = {name: fixture[name].copy() for name in fixture.files}

    for label in ("coarse", "fine"):
        assert source[f"{label}_phase"].shape == (1000,)
        assert source[f"{label}_wave"].shape == (1001, 2)
        for body in (1, 2):
            assert source[f"{label}_body{body}"].shape == (1001, 13)
        assert all(np.isfinite(source[f"{label}_{name}"]).all()
                   for name in ("phase", "wave", "body1", "body2"))

    np.testing.assert_array_equal(source["coarse_phase"], source["fine_phase"])
    np.testing.assert_array_equal(source["coarse_wave"], source["fine_wave"])
    time = np.arange(1001) * .01
    np.testing.assert_allclose(source["fine_wave"][:, 0], time, rtol=0, atol=1e-10)
    for body in (1, 2):
        coarse = source[f"coarse_body{body}"]
        fine = source[f"fine_body{body}"]
        np.testing.assert_allclose(coarse[:, 0], time, rtol=0, atol=1e-10)
        np.testing.assert_allclose(fine[:, 0], time, rtol=0, atol=1e-10)
        # This is a source-to-source negative control, not a parity pass.
        assert np.max(np.abs(coarse[:, 1] - fine[:, 1])) > .2

    phase_file = tmp_path / "phase.csv"
    np.savetxt(phase_file, source["fine_phase"], delimiter=",", fmt="%.17g")
    wec = WEC("Unmoored RM3 source-step diagnostic")
    float_body = wec.body("float", HYDRO, inertia=(0, 21_306_090.66, 0))
    spar = wec.body("spar", HYDRO, inertia=(0, 94_407_091.24, 0))
    wec.floating_joint(float_body, spar, damping=1_200_000)
    response = wec.run(
        JONSWAPWave(2, 8, phase_file=phase_file, discretization="traditional"),
        dt=.00125, end_time=10, ramp_time=0, radiation_memory=60,
        initial_coordinate={"spar_heave": -.21},
    )
    indices = np.arange(1001) * 8
    np.testing.assert_allclose(response.time[indices], time, rtol=0, atol=1e-9)
    np.testing.assert_allclose(response.wave_elevation[indices],
                               source["fine_wave"][:, 1], rtol=0, atol=1e-10)

    # Reuse the position and speed limits of the existing RM3 MoorDyn
    # MaxStep-only paired gate. No source motion or force enters this solve.
    limits = {"position": (0.002, 5e-5),
              "velocity": (0.0007, 5e-5)}
    for body, name in ((1, "float"), (2, "spar")):
        fine = source[f"fine_body{body}"]
        for kind, offset in (("position", 1), ("velocity", 7)):
            actual = getattr(response.bodies[name], kind)[indices]
            for axis in (0, 2, 4):
                error = float(np.max(np.abs(actual[:, axis] - fine[:, offset + axis])))
                assert error < limits[kind][axis == 4], (name, kind, axis, error)

"""Pack post-screen 150 s no-active-delay MaxStep controls by wave case.

Usage:
  python tools/pack_ellipsoid_ode45_continuous_nonlinear_mass_control_step.py \
      CASE ARTIFACT_ROOT OUTPUT.npz [--expect ARCHIVED.npz]

ARTIFACT_ROOT contains CASE_0p0125 and CASE_0p00625 directories downloaded
from the four named artifacts in GitHub Actions run 40131806849.
"""

import argparse
from pathlib import Path

import numpy as np

from pack_ellipsoid_ode45_continuous_nonlinear_mass_control import _check_hydro


ROOT = Path(__file__).resolve().parents[1]
BASE = (ROOT / "tests/fixtures/matlab"
        / "ellipsoid_ode45_continuous_nonlinear_mass_control")
STEPS = (("0p0125", .0125), ("0p00625", .00625))
ROWS = 3001
TIME = np.arange(ROWS) * .05


def _read(folder: Path, kind: str, shape: tuple[int, int]):
    values = np.loadtxt(folder / f"{kind}.csv", delimiter=",", ndmin=2)
    assert values.shape == shape and np.isfinite(values).all(), (folder, kind)
    return values


def pack(case: str, artifact_root: Path, output: Path,
         expected: Path | None = None) -> None:
    if case not in ("Regular", "RegularCIC"):
        raise ValueError(f"Unsupported ellipsoid wave case: {case}")
    basename = "regular_cic" if case == "RegularCIC" else "regular"
    with np.load(BASE / f"{basename}.npz", allow_pickle=False) as base:
        base_mass = base["mass"]
        base_wave = base["wave"]
        np.testing.assert_allclose(base["time"], TIME, rtol=0, atol=1e-9)
        np.testing.assert_allclose(base["settings"],
                                   [.05, .05, .05, .05, .05, 150, 0, 0, 0],
                                   rtol=0, atol=1e-12)
    packed: dict[str, np.ndarray] = {"case": np.array(case), "time": TIME}
    for tag, step in STEPS:
        artifact = artifact_root / f"{case}_{tag}"
        _check_hydro(artifact / "applications/Nonlinear_Hydro"
                     / "hydroData/ellipsoid.h5")
        folder = (artifact
                  / "matlab-ellipsoid-ode45-continuous-nonlinear-mass-control-step"
                  / case / f"maxstep_{tag}")
        body = _read(folder, "body", (ROWS, 55))
        pto = _read(folder, "pto", (ROWS, 25))
        wave = _read(folder, "wave", (ROWS, 2))
        applied = _read(folder, "applied_added", (ROWS, 7))
        settings = _read(folder, "settings", (1, 9))[0]
        mass = _read(folder, "mass", (1, 12))[0]
        np.testing.assert_allclose(settings,
                                   [.05, .05, .05, .05, step, 150, 0, 0, 0],
                                   rtol=0, atol=1e-12)
        for name, values in (("body", body), ("pto", pto), ("wave", wave),
                             ("applied", applied)):
            np.testing.assert_allclose(values[:, 0], TIME, rtol=0, atol=1e-9,
                                       err_msg=f"{case} {tag} {name} time")
        np.testing.assert_allclose(wave[:, 1], base_wave, rtol=0, atol=1e-10)
        np.testing.assert_allclose(np.delete(mass, 6),
                                   np.delete(base_mass, 6), rtol=0, atol=1e-7)
        assert abs(mass[6] - step) < 1e-12
        inactive_speed = body[:, [7, 8, 10, 11, 12]]
        inactive_acceleration = body[:, [49, 50, 52, 53, 54]]
        assert np.max(np.abs(inactive_speed)) < 1e-9
        assert np.max(np.abs(inactive_acceleration)) < 1e-9
        assert np.max(np.abs(inactive_acceleration @ mass[7:])) < 1e-6
        assert np.max(np.abs(applied[:, 3])) < 1e-6
        np.testing.assert_allclose(pto[:, 3], body[:, 3] + 2,
                                   rtol=0, atol=1e-9)
        np.testing.assert_allclose(pto[:, 9], body[:, 9], rtol=0, atol=1e-9)
        np.testing.assert_allclose(pto[:, 15], -1_200_000 * pto[:, 9],
                                   rtol=0, atol=1e-6)
        np.testing.assert_allclose(pto[:, 21], pto[:, 15] * pto[:, 9],
                                   rtol=0, atol=1e-6)
        prefix = f"maxstep_{tag}_"
        packed.update({
            prefix + "heave": body[:, 3],
            prefix + "speed": body[:, 9],
            prefix + "total": body[:, 15],
            prefix + "excitation": body[:, 21],
            prefix + "radiation": body[:, 27],
            prefix + "reported_added_mass": body[:, 33],
            prefix + "restoring": body[:, 39],
            prefix + "viscous": body[:, 45],
            prefix + "acceleration": body[:, 51],
            prefix + "applied_heave_added": applied[:, 3],
            prefix + "inactive_speed": inactive_speed,
            prefix + "inactive_acceleration": inactive_acceleration,
            prefix + "pto_stroke": pto[:, 3],
            prefix + "pto_speed": pto[:, 9],
            prefix + "pto_force": pto[:, 15],
            prefix + "pto_power": -pto[:, 21],
            prefix + "wave": wave[:, 1],
            prefix + "settings": settings,
            prefix + "mass": mass,
        })
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output, **packed)
    print(f"Packed {case} no-active-delay MaxStep traces to {output}")
    if expected is not None:
        with np.load(expected, allow_pickle=False) as archived:
            assert set(archived.files) == set(packed), expected
            for name, value in packed.items():
                np.testing.assert_equal(value, archived[name],
                                        err_msg=f"{case} {name} changed")
        print(f"Matched archived no-active-delay traces in {expected}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", choices=("Regular", "RegularCIC"))
    parser.add_argument("artifact_root", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--expect", type=Path)
    args = parser.parse_args()
    pack(args.case, args.artifact_root, args.output, args.expect)

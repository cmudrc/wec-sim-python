"""Pack copied-model ellipsoid MaxStep artifacts into a portable fixture.

Usage:
  python tools/pack_ellipsoid_ode45_restoring_step.py CASE ARTIFACT_ROOT OUT.npz
      [--expect ARCHIVED.npz]

ARTIFACT_ROOT contains the eight named artifacts from run 40131187182.
The archived 0.05 s continuous-nonlinear arm supplies the common base trace.
"""

import argparse
from pathlib import Path

import numpy as np

from pack_ellipsoid_ode45_restoring_timing import _check_source_hydro


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "tests/fixtures/matlab/ellipsoid_ode45_restoring_timing"
STEPS = (("0.025", "0p025", .025),
         ("0.0125", "0p0125", .0125),
         ("0.00625", "0p00625", .00625),
         ("0.003125", "0p003125", .003125))
ROWS = 3001
TIME = np.arange(ROWS) * .05


def _read(folder: Path, kind: str, shape: tuple[int, int]) -> np.ndarray:
    values = np.loadtxt(folder / f"continuous_nonlinear_{kind}.csv",
                        delimiter=",", ndmin=2)
    assert values.shape == shape and np.isfinite(values).all(), (folder, kind)
    return values


def pack(case: str, artifact_root: Path, output: Path,
         expected: Path | None = None) -> None:
    if case not in ("Regular", "RegularCIC"):
        raise ValueError(f"Unsupported ellipsoid wave case: {case}")
    basename = "regular_cic" if case == "RegularCIC" else "regular"
    with np.load(BASE / f"{basename}.npz", allow_pickle=False) as base:
        assert str(base["case"]) == case
        base_time = base["continuous_nonlinear_time"]
        base_wave = base["continuous_nonlinear_wave"]
        np.testing.assert_allclose(
            base["continuous_nonlinear_settings"],
            [.05, .05, .05, .05, .05, 150, 0, 0, 0], rtol=0, atol=1e-12,
        )
    np.testing.assert_allclose(base_time, TIME, rtol=0, atol=1e-9)
    packed: dict[str, np.ndarray] = {
        "case": np.array(case),
        "time": base_time,
    }
    for step_text, tag, step in STEPS:
        artifact = (artifact_root / f"ellipsoid-ode45-restoring-{case}"
                    f"-maxstep-{step_text}")
        _check_source_hydro(
            artifact / "applications/Nonlinear_Hydro/hydroData/ellipsoid.h5"
        )
        folder = (artifact / "matlab-ellipsoid-ode45-restoring-timing-step"
                  / case / f"maxstep_{tag}")
        body = _read(folder, "body", (ROWS, 55))
        pto = _read(folder, "pto", (ROWS, 25))
        wave = _read(folder, "wave", (ROWS, 2))
        settings = _read(folder, "settings", (1, 9))[0]
        np.testing.assert_allclose(settings,
                                   [.05, .05, .05, .05, step, 150, 0, 0, 0],
                                   rtol=0, atol=1e-12,
                                   err_msg=f"{case} {step_text} settings")
        for name, values in (("body", body), ("pto", pto), ("wave", wave)):
            np.testing.assert_allclose(values[:, 0], TIME, rtol=0, atol=1e-9,
                                       err_msg=f"{case} {step_text} {name} time")
        np.testing.assert_allclose(wave[:, 1], base_wave, rtol=0, atol=1e-10)
        np.testing.assert_allclose(pto[:, 3], body[:, 3] + 2,
                                   rtol=0, atol=1e-10)
        np.testing.assert_allclose(pto[:, 9], body[:, 9], rtol=0, atol=1e-9)
        np.testing.assert_allclose(pto[:, 15], -1_200_000 * pto[:, 9],
                                   rtol=0, atol=1e-6)
        np.testing.assert_allclose(pto[:, 21], pto[:, 15] * pto[:, 9],
                                   rtol=0, atol=1e-6)
        prefix = f"maxstep_{tag}_"
        packed.update({
            prefix + "heave": body[:, 3],
            prefix + "speed": body[:, 9],
            prefix + "restoring": body[:, 39],
            prefix + "pto_stroke": pto[:, 3],
            prefix + "pto_speed": pto[:, 9],
            prefix + "pto_force": pto[:, 15],
            prefix + "pto_power": -pto[:, 21],
            prefix + "wave": wave[:, 1],
            prefix + "settings": settings,
        })
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output, **packed)
    print(f"Packed {case} MaxStep traces to {output}")
    if expected is not None:
        with np.load(expected, allow_pickle=False) as archived:
            assert set(archived.files) == set(packed), expected
            for name, value in packed.items():
                np.testing.assert_equal(
                    value, archived[name],
                    err_msg=f"{case} {name} differs from {expected}",
                )
        print(f"Matched archived MaxStep traces in {expected}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", choices=("Regular", "RegularCIC"))
    parser.add_argument("artifact_root", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--expect", type=Path)
    args = parser.parse_args()
    pack(args.case, args.artifact_root, args.output, args.expect)

"""Pack one pinned source job's corrected ellipsoid traces for portable CI.

Usage:
  python tools/pack_ellipsoid_ode45_mass_correction.py CASE CSV_DIR OUTPUT.npz
  [--expect ARCHIVED.npz]

CASE is Regular or RegularCIC. CSV_DIR is the case folder uploaded by the
Ellipsoid ode45 derived heave mass correction workflow.
"""

import argparse
from pathlib import Path

import numpy as np


SOURCE_SHAPES = {"body": (3001, 55), "pto": (3001, 25),
                 "wave": (3001, 2), "applied_added": (3001, 7),
                 "mass": (1, 12)}


def _read(folder: Path, mode: str, name: str) -> np.ndarray:
    values = np.loadtxt(folder / f"{mode}_{name}.csv", delimiter=",", ndmin=2)
    assert values.shape == SOURCE_SHAPES[name], (mode, name, values.shape)
    assert np.isfinite(values).all(), (mode, name)
    return values


def pack(case: str, folder: Path, output: Path,
         expected: Path | None = None) -> None:
    if case not in {"Regular", "RegularCIC"}:
        raise ValueError(f"Unsupported wave case: {case}")
    packed = {"case": np.array(case)}
    common_time = np.arange(3001) * .05
    for mode, expected_step in (("published_clock", .05),
                                ("refined_clock", .005)):
        body, pto, wave, applied, mass = (
            _read(folder, mode, name) for name in SOURCE_SHAPES
        )
        for name, record in (("body", body), ("pto", pto),
                             ("wave", wave), ("applied_added", applied)):
            np.testing.assert_allclose(record[:, 0], common_time,
                                       rtol=0, atol=1e-10, err_msg=name)
        assert abs(mass[0, 6] - expected_step) < 1e-12
        prefix = f"{mode}_"
        packed.update({
            prefix + "time": common_time,
            prefix + "heave": body[:, 3],
            prefix + "speed": body[:, 9],
            prefix + "pto_force": pto[:, 15],
            prefix + "pto_speed": pto[:, 9],
            prefix + "wave": wave[:, 1],
            prefix + "applied_heave": applied[:, 3],
            prefix + "inactive_speed": body[:, [7, 8, 10, 11, 12]],
            prefix + "inactive_acceleration": body[:, [49, 50, 52, 53, 54]],
            prefix + "mass": mass[0],
        })
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output, **packed)
    print(f"Packed {case} corrected source to {output}")
    if expected is not None:
        with np.load(expected, allow_pickle=False) as archive:
            assert set(archive.files) == set(packed), expected
            for name, values in packed.items():
                if name == "case":
                    assert str(archive[name]) == str(values)
                else:
                    np.testing.assert_allclose(
                        values, archive[name], rtol=1e-9, atol=1e-9,
                        err_msg=f"{case} {name} differs from {expected}",
                    )
        print(f"Matched archived source traces in {expected}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", choices=("Regular", "RegularCIC"))
    parser.add_argument("csv_dir", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--expect", type=Path)
    args = parser.parse_args()
    pack(args.case, args.csv_dir, args.output, args.expect)

"""Pack a full 150 s ellipsoid continuous-nonlinear mass-control artifact.

Usage:
  python tools/pack_ellipsoid_ode45_continuous_nonlinear_mass_control.py \
      CASE CSV_DIR OUTPUT.npz [--expect ARCHIVED.npz] [--source-hydro H5]

The MATLAB workflow copies the application input and model. This script reads
its output, checks the matched clocks and no-active-heave-feedback identity,
and stores the source traces needed by the portable diagnostic gate.
"""

import argparse
from pathlib import Path

import h5py
import numpy as np


ROWS = 3001
TIMES = np.arange(ROWS) * .05
ROOT = Path(__file__).resolve().parents[1]
PINNED_HYDRO = (ROOT / "tests/fixtures/matlab/ellipsoid_ode45_convergence"
                / "hydroData/ellipsoid.h5")


def _read(folder: Path, name: str, shape: tuple[int, int]) -> np.ndarray:
    values = np.loadtxt(folder / f"{name}.csv", delimiter=",", ndmin=2)
    assert values.shape == shape, (name, values.shape)
    assert np.isfinite(values).all(), name
    return values


def _datasets(group, prefix=""):
    result = {}
    for name, item in group.items():
        path = f"{prefix}/{name}"
        if isinstance(item, h5py.Dataset):
            result[path] = item[()]
        else:
            result.update(_datasets(item, path))
    return result


def _check_hydro(path: Path) -> None:
    with h5py.File(path) as fresh, h5py.File(PINNED_HYDRO) as archived:
        fresh_data, archived_data = _datasets(fresh), _datasets(archived)
    assert len(fresh_data) == 38 and fresh_data.keys() == archived_data.keys()
    for name in fresh_data:
        np.testing.assert_equal(fresh_data[name], archived_data[name],
                                err_msg=f"source HDF5 dataset {name} changed")
    print(f"Matched 38 pinned HDF5 datasets in {path}")


def pack(case: str, folder: Path, output: Path,
         expected: Path | None = None, source_hydro: Path | None = None) -> None:
    if case not in ("Regular", "RegularCIC"):
        raise ValueError(f"Unsupported ellipsoid wave case: {case}")
    if source_hydro is not None:
        _check_hydro(source_hydro)

    body = _read(folder, "body", (ROWS, 55))
    pto = _read(folder, "pto", (ROWS, 25))
    wave = _read(folder, "wave", (ROWS, 2))
    applied = _read(folder, "applied_added", (ROWS, 7))
    settings = _read(folder, "settings", (1, 9))[0]
    mass = _read(folder, "mass", (1, 12))[0]
    np.testing.assert_allclose(settings, [.05]*5 + [150, 0, 0, 0],
                               rtol=0, atol=1e-12)
    for name, values in (("body", body), ("pto", pto), ("wave", wave),
                         ("applied", applied)):
        np.testing.assert_allclose(values[:, 0], TIMES, rtol=0, atol=1e-9,
                                   err_msg=f"{name} changed its output grid")
    np.testing.assert_allclose(mass[0], mass[2] / mass[1], rtol=0, atol=1e-13)
    np.testing.assert_allclose(mass[4], mass[3] + mass[2], rtol=0, atol=1e-7)
    assert abs(mass[5]) < 1e-7 and abs(mass[6] - .05) < 1e-12
    inactive_speed = body[:, [7, 8, 10, 11, 12]]
    inactive_acceleration = body[:, [49, 50, 52, 53, 54]]
    assert np.max(np.abs(inactive_speed)) < 1e-9
    assert np.max(np.abs(inactive_acceleration)) < 1e-9
    assert np.max(np.abs(inactive_acceleration @ mass[7:])) < 1e-6
    assert np.max(np.abs(applied[:, 3])) < 1e-6
    np.testing.assert_allclose(pto[:, 3], body[:, 3] + 2, rtol=0, atol=1e-9)
    np.testing.assert_allclose(pto[:, 9], body[:, 9], rtol=0, atol=1e-9)
    np.testing.assert_allclose(pto[:, 15], -1_200_000 * pto[:, 9],
                               rtol=0, atol=1e-6)
    np.testing.assert_allclose(pto[:, 21], pto[:, 15] * pto[:, 9],
                               rtol=0, atol=1e-6)

    packed = {
        "case": np.array(case), "time": body[:, 0],
        "heave": body[:, 3], "speed": body[:, 9],
        "total": body[:, 15], "excitation": body[:, 21],
        "radiation": body[:, 27], "reported_added_mass": body[:, 33],
        "restoring": body[:, 39], "viscous": body[:, 45],
        "acceleration": body[:, 51],
        "applied_heave_added": applied[:, 3],
        "inactive_speed": inactive_speed,
        "inactive_acceleration": inactive_acceleration,
        "pto_stroke": pto[:, 3], "pto_speed": pto[:, 9],
        "pto_force": pto[:, 15], "pto_power": -pto[:, 21],
        "wave": wave[:, 1], "settings": settings, "mass": mass,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output, **packed)
    print(f"Packed {case} full-duration mass-control traces to {output}")
    if expected is not None:
        with np.load(expected, allow_pickle=False) as archived:
            assert set(archived.files) == set(packed)
            for name, value in packed.items():
                if name == "case":
                    assert str(archived[name]) == case
                else:
                    np.testing.assert_allclose(
                        value, archived[name], rtol=1e-9, atol=1e-9,
                        err_msg=f"{case} {name} differs from {expected}",
                    )
        print(f"Matched archived mass-control traces in {expected}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", choices=("Regular", "RegularCIC"))
    parser.add_argument("csv_dir", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--expect", type=Path)
    parser.add_argument("--source-hydro", type=Path)
    args = parser.parse_args()
    pack(args.case, args.csv_dir, args.output, args.expect, args.source_hydro)

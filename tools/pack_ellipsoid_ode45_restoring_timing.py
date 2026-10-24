"""Pack one full-duration ellipsoid force-timing source artifact for Python CI.

Usage:
  python tools/pack_ellipsoid_ode45_restoring_timing.py CASE CSV_DIR OUTPUT.npz
      [--expect ARCHIVED.npz] [--source-hydro FRESH.h5]

CSV_DIR contains the three 150 s arms uploaded by the repository-owned
ellipsoid-ode45-restoring-timing workflow. The original MATLAB repositories
are never modified by this tool.
"""

import argparse
from pathlib import Path

import h5py
import numpy as np


MODES = ("baseline", "continuous_restoring", "continuous_nonlinear")
ROWS = 3001
DT = .05
EXPECTED_RATES = {
    "baseline": (.05, .05, .05),
    "continuous_restoring": (0., .05, .05),
    "continuous_nonlinear": (0., 0., 0.),
}


def _read(folder: Path, mode: str, name: str, shape: tuple[int, int]):
    values = np.loadtxt(folder / f"{mode}_{name}.csv", delimiter=",", ndmin=2)
    assert values.shape == shape, (mode, name, values.shape)
    assert np.isfinite(values).all(), (mode, name)
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


def _check_source_hydro(path: Path) -> None:
    reference = (Path(__file__).resolve().parents[1]
                 / "tests/fixtures/matlab/ellipsoid_ode45_convergence"
                 / "hydroData/ellipsoid.h5")
    with h5py.File(path) as source, h5py.File(reference) as archived:
        source_data, archived_data = _datasets(source), _datasets(archived)
    assert len(source_data) == 38 and source_data.keys() == archived_data.keys()
    for name in source_data:
        np.testing.assert_equal(source_data[name], archived_data[name],
                                err_msg=f"source HDF5 dataset {name} changed")
    print(f"Matched 38 pinned HDF5 datasets in {path}")


def pack(case: str, folder: Path, output: Path,
         expected: Path | None = None,
         source_hydro: Path | None = None) -> None:
    if case not in ("Regular", "RegularCIC"):
        raise ValueError(f"Unsupported ellipsoid wave case: {case}")
    if source_hydro is not None:
        _check_source_hydro(source_hydro)
    times = np.arange(ROWS) * DT
    packed: dict[str, np.ndarray] = {"case": np.array(case)}
    waves = []
    for mode in MODES:
        body = _read(folder, mode, "body", (ROWS, 55))
        pto = _read(folder, mode, "pto", (ROWS, 25))
        wave = _read(folder, mode, "wave", (ROWS, 2))
        settings = _read(folder, mode, "settings", (1, 9))[0]
        np.testing.assert_allclose(
            settings, [.05, .05, .05, .05, .05, 150, *EXPECTED_RATES[mode]],
            rtol=0, atol=1e-12, err_msg=f"{mode} changed a pinned clock",
        )
        for name, source in (("body", body), ("pto", pto), ("wave", wave)):
            np.testing.assert_allclose(source[:, 0], times, rtol=0,
                                       atol=1e-9, err_msg=f"{mode} {name} time")
        np.testing.assert_allclose(pto[:, 9], body[:, 9], rtol=0, atol=1e-9,
                                   err_msg=f"{mode} inactive PTO speed")
        np.testing.assert_allclose(pto[:, 15], -1_200_000 * pto[:, 9],
                                   rtol=0, atol=1e-6,
                                   err_msg=f"{mode} PTO damping law")
        np.testing.assert_allclose(pto[:, 21], pto[:, 15] * pto[:, 9],
                                   rtol=0, atol=1e-6,
                                   err_msg=f"{mode} PTO source power law")
        waves.append(wave[:, 1])
        prefix = f"{mode}_"
        packed.update({
            prefix + "time": body[:, 0],
            prefix + "heave": body[:, 3],
            prefix + "speed": body[:, 9],
            prefix + "restoring": body[:, 39],
            prefix + "excitation": body[:, 21],
            prefix + "acceleration": body[:, 51],
            prefix + "pto_stroke": pto[:, 3],
            prefix + "pto_speed": pto[:, 9],
            prefix + "pto_force": pto[:, 15],
            prefix + "pto_power": -pto[:, 21],
            prefix + "wave": wave[:, 1],
            prefix + "settings": settings,
        })
    for wave in waves[1:]:
        np.testing.assert_allclose(wave, waves[0], rtol=0, atol=1e-10,
                                   err_msg="force timing changed the incident wave")
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output, **packed)
    print(f"Packed {case} source timing traces to {output}")
    if expected is not None:
        with np.load(expected, allow_pickle=False) as archived:
            assert set(archived.files) == set(packed), expected
            for name, value in packed.items():
                if name == "case":
                    assert str(archived[name]) == case
                else:
                    np.testing.assert_allclose(
                        value, archived[name], rtol=1e-9, atol=1e-9,
                        err_msg=f"{case} {name} differs from {expected}",
                    )
        print(f"Matched archived timing traces in {expected}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", choices=("Regular", "RegularCIC"))
    parser.add_argument("csv_dir", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--expect", type=Path)
    parser.add_argument("--source-hydro", type=Path)
    args = parser.parse_args()
    pack(args.case, args.csv_dir, args.output, args.expect,
         args.source_hydro)

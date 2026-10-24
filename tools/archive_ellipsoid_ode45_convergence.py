"""Archive full 150 s ellipsoid ode45 source-screen artifacts for portable CI.

Example:
  gh run download RUN_ID --dir /tmp/ellipsoid-source
  python tools/archive_ellipsoid_ode45_convergence.py \
      /tmp/ellipsoid-source tests/fixtures/matlab/ellipsoid_ode45_convergence
"""

import argparse
import hashlib
from pathlib import Path
import shutil

import h5py
import numpy as np
from scipy.integrate import trapezoid


ROOT = Path(__file__).resolve().parents[1]
ARCHIVED_005 = (ROOT / "tests/fixtures/matlab/ellipsoid_ode45_maxstep_only"
                / "refined_005_common_grid.npz")
STEPS = {"step005": .005, "step0025": .0025,
         "step00125": .00125, "step000625": .000625}
OUTPUT_STEPS = {label: (.005 if label == "step000625" else step)
                for label, step in STEPS.items()}
STABILITY_LIMITS = {
    "Regular": (.001, .0013, 1600, 1300),
    "RegularCIC": (.0011, .0015, 2000, 1600),
}
CHANNELS = (1, 2, 5, 6)
WORK_STEP_LIMIT = .001
OUTPUT_WORK_LIMIT = 1e-5
IRF_DATASETS = {
    "body1/hydro_coeffs/excitation/impulse_response_fun/f",
    "body1/hydro_coeffs/radiation_damping/impulse_response_fun/K",
}


def _sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _maximum(a, b):
    assert a.shape == b.shape
    assert np.isfinite(a).all() and np.isfinite(b).all()
    return float(np.max(np.abs(a - b)))


def _datasets(path):
    result = {}
    with h5py.File(path) as source:
        def collect(name, obj):
            if isinstance(obj, h5py.Dataset):
                result[name] = obj[()]
        source.visititems(collect)
    return result


def _read_case(artifacts, case):
    base = artifacts / f"ellipsoid-ode45-convergence-{case}"
    output = base / "matlab-ellipsoid-ode45-convergence" / case
    histories = {}
    settings_records = {}
    for label, step in STEPS.items():
        output_step = OUTPUT_STEPS[label]
        history = np.loadtxt(output / f"{label}_history.csv",
                             delimiter=",", ndmin=2)
        settings = np.loadtxt(output / f"{label}_settings.csv",
                              delimiter=",", ndmin=2)
        assert history.shape == (round(150 / output_step) + 1, 11)
        assert np.isfinite(history).all()
        np.testing.assert_allclose(settings,
                                   [[step, output_step, step, step, step]],
                                   rtol=0, atol=1e-12)
        np.testing.assert_allclose(history[:, 0],
                                   np.arange(len(history)) * output_step,
                                   rtol=0, atol=1e-9)
        # The PTO is constrained to vertical heave. This identity guards
        # against accidentally archiving the inactive local PTO axis.
        assert _maximum(history[:, 4], history[:, 2]) < 1e-9
        assert _maximum(history[:, 5], -1_200_000 * history[:, 4]) < 1e-6
        assert _maximum(history[:, 6], -history[:, 5] * history[:, 4]) < 1e-6
        histories[label] = history
        settings_records[label] = settings
    return base, histories, settings_records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifacts", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    records = {}
    inputs = {}
    with np.load(ARCHIVED_005, allow_pickle=False) as archived:
        for case in ("Regular", "RegularCIC"):
            base, histories, settings_records = _read_case(args.artifacts, case)
            key = "regular_cic" if case == "RegularCIC" else "regular"
            five = histories["step005"]
            for column, suffix, limit in ((0, "time", 1e-10),
                                          (1, "heave", 1e-7),
                                          (2, "speed", 1e-7),
                                          (5, "pto_force", .1),
                                          (7, "wave", 1e-10),
                                          (8, "restoring", .1)):
                assert _maximum(five[::10, column],
                                archived[f"{key}_{suffix}"]) < limit
            half = histories["step0025"]
            quarter = histories["step00125"]
            finest = histories["step000625"]
            assert _maximum(five[:, 7], half[::2, 7]) < 1e-10
            assert _maximum(half[:, 7], quarter[::2, 7]) < 1e-10
            assert _maximum(quarter[::4, 7], finest[:, 7]) < 1e-10
            common = (five, half[::2], quarter[::4], finest)
            for column, limit in zip(CHANNELS, STABILITY_LIMITS[case],
                                     strict=True):
                first = _maximum(common[0][:, column], common[1][:, column])
                previous = _maximum(common[1][:, column], common[2][:, column])
                final = _maximum(common[2][:, column], common[3][:, column])
                print(case, column, "source changes", first, previous, final)
                assert (final < limit and previous < .7 * first
                        and final < .7 * previous)
            works = [float(trapezoid(history[:, 6], history[:, 0]))
                     for history in (five, half, quarter, finest)]
            assert all(np.isfinite(work) and work > 1e3 for work in works)
            common_works = [float(trapezoid(history[:, 6], history[:, 0]))
                            for history in common]
            for label, full, sampled in zip(
                list(STEPS)[:3], works[:3], common_works[:3], strict=True
            ):
                error = abs(full / sampled - 1)
                print(case, label, "full/common signed-work difference", error)
                assert error < OUTPUT_WORK_LIMIT
            work_changes = [abs(later - earlier) / later
                            for earlier, later in zip(works, works[1:])]
            common_work_changes = [abs(later - earlier) / later
                                   for earlier, later in zip(
                                       common_works, common_works[1:])]
            print(case, "source signed PTO work (J)", works)
            print(case, "source signed PTO work relative changes", work_changes)
            print(case, "common-grid signed PTO work (J)", common_works)
            print(case, "common-grid work relative changes",
                  common_work_changes)
            for changes in (work_changes, common_work_changes):
                assert (changes[-1] < WORK_STEP_LIMIT
                        and changes[1] < .7 * changes[0]
                        and changes[2] < .7 * changes[1])
            for label, history in histories.items():
                records[f"{key}_{label}"] = history
                records[f"{key}_{label}_settings"] = settings_records[label]
            inputs[case] = {
                "hydro": base / "applications/Nonlinear_Hydro/hydroData/ellipsoid.h5",
                "geometry": base / "applications/Nonlinear_Hydro/geometry/elipsoid.stl",
            }

    assert _sha(inputs["Regular"]["geometry"]) == _sha(
        inputs["RegularCIC"]["geometry"])
    regular_hydro = _datasets(inputs["Regular"]["hydro"])
    cic_hydro = _datasets(inputs["RegularCIC"]["hydro"])
    assert regular_hydro.keys() == cic_hydro.keys()
    for name, values in regular_hydro.items():
        other = cic_hydro[name]
        if name in IRF_DATASETS:
            # Independent bemio runs differ only at roundoff in the two
            # computed impulse-response arrays (observed <2.1e-12 absolute
            # in the prior pinned R2025b artifacts). Keep all other BEM
            # datasets byte-for-byte equal in value.
            np.testing.assert_allclose(values, other, rtol=1e-12,
                                       atol=1e-11, equal_nan=True,
                                       err_msg=name)
        else:
            np.testing.assert_array_equal(values, other, err_msg=name)

    args.output.mkdir(parents=True, exist_ok=True)
    destination = args.output / "source.npz"
    np.savez_compressed(destination, **records)
    (args.output / "hydroData").mkdir(exist_ok=True)
    (args.output / "geometry").mkdir(exist_ok=True)
    shutil.copyfile(inputs["Regular"]["hydro"],
                    args.output / "hydroData/ellipsoid.h5")
    shutil.copyfile(inputs["Regular"]["geometry"],
                    args.output / "geometry/elipsoid.stl")
    for file in (destination, args.output / "hydroData/ellipsoid.h5",
                 args.output / "geometry/elipsoid.stl"):
        print(file, _sha(file), file.stat().st_size)


if __name__ == "__main__":
    main()

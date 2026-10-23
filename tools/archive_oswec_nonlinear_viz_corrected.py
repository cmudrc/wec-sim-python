"""Archive the pinned corrected OSWEC visualization run for ordinary CI."""

from argparse import ArgumentParser
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile, ZipInfo

import numpy as np
from scipy.io import loadmat


def archive(source_dir: Path, output_dir: Path) -> None:
    baseline = loadmat(source_dir / "baseline.mat")
    stability = loadmat(source_dir / "stability.mat")
    assert float(baseline["dt"].item()) == 0.0125
    assert float(stability["dt"].item()) == 0.00625
    time = baseline["time"].ravel()
    np.testing.assert_allclose(time, np.arange(1201) * 0.1, rtol=0, atol=1e-11)
    np.testing.assert_allclose(stability["time"].ravel(), time, rtol=0, atol=1e-11)
    for key in ("ptoTime", "pressureTime"):
        np.testing.assert_allclose(baseline[key].ravel(), time, rtol=0, atol=1e-11)
    np.testing.assert_allclose(stability["ptoTime"].ravel(), time,
                               rtol=0, atol=1e-11)

    arrays = {"time": time, "source_dt": baseline["dt"].ravel(),
              "stability_dt": stability["dt"].ravel(),
              "wave_number": baseline["waveNumber"].ravel(),
              "water_depth": baseline["waterDepth"].ravel(),
              "areas": baseline["areas"].ravel()}
    for key in ("pose", "velocity", "basePose", "baseVelocity",
                "ptoPosition", "ptoVelocity", "ptoForce", "ptoPower"):
        value = baseline[key]
        assert value.shape == (1201, 6) and np.isfinite(value).all(), key
        arrays[key] = value
        arrays[f"stable_{key}"] = stability[key]
        assert stability[key].shape == value.shape
    for key in ("hydrostatic", "nonlinearWave", "linearWave"):
        value = baseline[key]
        assert value.shape == (1201, 1042) and np.isfinite(value).all(), key
        arrays[key] = value
    assert arrays["areas"].shape == (1042,)
    output_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output_dir / "history.npz", **arrays)

    samples = source_dir / "vtp_samples"
    names = ("ground.txt", "bodies.txt", "OSWEC.pvd", *(
        f"{folder}/{stem}_{frame}.vtp"
        for folder, stem in (("waves", "waves"),
                             ("body1_flap", "flap"),
                             ("body2_base", "base"))
        for frame in (1, 601, 1201)
    ))
    with ZipFile(output_dir / "vtp_samples.zip", "w", ZIP_DEFLATED, 9) as archive_file:
        for name in names:
            info = ZipInfo(name, date_time=(2026, 1, 1, 0, 0, 0))
            info.compress_type = ZIP_DEFLATED
            archive_file.writestr(info, (samples / name).read_bytes())


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("source_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    archive(args.source_dir, args.output_dir)

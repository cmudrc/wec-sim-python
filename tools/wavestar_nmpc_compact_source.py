"""Pack one checked WaveStar source artifact for portable diagnostics."""

import argparse
import json
from pathlib import Path

import numpy as np

from wavestar_nmpc_reference_compare import (
    CASE, FIELDS, PREFIX, ROOT, _read_artifact,
)

EXPECTED_SOURCE = {
    "wec_sim_commit": "0753b2e47f2457c078751dcfe5d251d1767b80ab",
    "applications_commit": "d53d4d4c9eda2581f04204f5d394a6ef84bb099e",
    "hydro_sha256": "b2541ff01c319987d9157d1be9fbc9ec39a7660d0572d8e316c428b8aa3359d1",
    "plant_step_s": .00025,
    "controller_step_s": .05,
    "samples": 4501,
}


def pack(artifact: Path, output: Path, step: float = .00025):
    if step not in (.00025, .000125):
        raise ValueError("Unsupported WaveStar source plant step")
    trace, _ = _read_artifact(artifact, step)
    audit = json.loads((artifact / "wavestar-reference-audit.json").read_text())
    expected_source = {**EXPECTED_SOURCE, "plant_step_s": step}
    if audit["source"] != expected_source:
        raise ValueError("Source provenance does not match the pinned case")
    components = np.loadtxt(artifact / f"{PREFIX}components.csv",
                            delimiter=",", ndmin=2)
    pto = np.loadtxt(artifact / f"{PREFIX}{CASE}_pto1.csv",
                     delimiter=",", ndmin=2)
    if components.shape != (500, 4) or pto.shape != (4501, 43):
        raise ValueError("Unexpected source components or PTO shape")
    if not np.isfinite(components).all() or not np.isfinite(pto).all():
        raise ValueError("Nonfinite source components or PTO values")
    with np.load(ROOT / "tests/fixtures/matlab/"
                 "wavestar_nmpc_split_step_225s/source.npz",
                 allow_pickle=False) as pinned:
        np.testing.assert_allclose(components, pinned["components"],
                                   rtol=0, atol=1e-12)
    np.savez_compressed(
        output,
        time=np.arange(4501) * .05,
        components=components,
        pto_speed=pto[:, 9],
        **{key: trace[key] for key in FIELDS.values()},
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--step", type=float, default=.00025)
    args = parser.parse_args()
    pack(args.artifact, args.output, args.step)


if __name__ == "__main__":
    main()

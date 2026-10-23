"""Archive verified MATLAB traces from one continuous-yaw CI run."""

import argparse
import hashlib
import json
import shutil
from pathlib import Path

from oswec_continuous_yaw_step_audit import audit


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def unique_file(root, name):
    matches = list(root.rglob(name))
    if len(matches) != 1:
        raise ValueError(f"Expected one {name} below {root}; found {len(matches)}")
    return matches[0]


def archive(input_root, output_dir, run_id, head_sha, patch):
    if output_dir.exists():
        raise FileExistsError(output_dir)
    manifest = {
        "run_url": f"https://github.com/cmudrc/wec-sim-python/actions/runs/{run_id}",
        "head_sha": head_sha,
        "continuous_heading_patch_sha256": sha256(patch),
        "seeds": {},
    }
    files = []
    for seed in (1, 2, 3):
        summary_path = unique_file(input_root, f"seed-{seed}-summary.json")
        source_dir = summary_path.parent
        report = json.loads(summary_path.read_text())
        calculated = audit(source_dir, seed)
        if report != calculated or not calculated["within_limits"]:
            raise ValueError(f"Seed {seed} source audit mismatch or failed limits")
        manifest["seeds"][str(seed)] = calculated
        for dt in (0.01, 0.005):
            source = source_dir / f"seed-{seed}-dt-{dt:g}.mat"
            if sha256(source) != calculated["source_sha256"][f"{dt:g}"]:
                raise ValueError(f"Seed {seed}, dt {dt:g} hash mismatch")
            files.append(source)
    output_dir.mkdir(parents=True)
    for source in files:
        shutil.copyfile(source, output_dir / source.name)
    (output_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--run-id", required=True, type=int)
    parser.add_argument("--head-sha", required=True)
    parser.add_argument("--patch", required=True, type=Path)
    args = parser.parse_args()
    archive(args.input_root, args.output_dir, args.run_id,
            args.head_sha, args.patch)


if __name__ == "__main__":
    main()

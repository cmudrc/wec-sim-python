"""Compare two independently generated WaveStar NMPC source references."""

import argparse
import json
from pathlib import Path

from wavestar_nmpc_reference_compare import _compare, _read_artifact


def compare(left_dir, right_dir):
    left_audit = json.loads((left_dir / "wavestar-reference-audit.json").read_text())
    right_audit = json.loads((right_dir / "wavestar-reference-audit.json").read_text())
    if left_audit["source"] != right_audit["source"]:
        raise ValueError("The repeated runs did not share source inputs and limits")
    step = left_audit["source"]["plant_step_s"]
    if step not in (.00025, .000125, .0000625):
        raise ValueError("Unsupported WaveStar source plant step")
    left, left_limits = _read_artifact(left_dir, step)
    right, right_limits = _read_artifact(right_dir, step)
    if left_limits != right_limits:
        raise ValueError("The repeated runs did not share pointwise limits")
    result = {
        "source": left_audit["source"],
        "pointwise_limits": left_limits,
        "repeat_comparison": _compare(left, right, left_limits),
        "python_physical_pointwise_left": left_audit["python_physical_pointwise"],
        "python_physical_pointwise_right": right_audit["python_physical_pointwise"],
    }
    for label, directory in (("left", left_dir), ("right", right_dir)):
        runtime = directory / "matlab-runtime.txt"
        cpu = directory / "cpu.txt"
        result[f"matlab_runtime_{label}"] = (
            runtime.read_text() if runtime.exists() else None
        )
        result[f"cpu_{label}"] = cpu.read_text() if cpu.exists() else None
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--left", type=Path, required=True)
    parser.add_argument("--right", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = compare(args.left, args.right)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({
        "source": result["source"],
        "repeat_comparison": result["repeat_comparison"],
        "python_physical_pointwise_left": result["python_physical_pointwise_left"],
        "python_physical_pointwise_right": result["python_physical_pointwise_right"],
    }, indent=2))


if __name__ == "__main__":
    main()

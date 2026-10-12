"""Render the published-case motion badge from the pinned evidence ledger."""

import argparse
import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
INVENTORY = ROOT / "tests/reference_cases.csv"
STATUS = ROOT / "tests/reference_case_status.csv"
BADGE = ROOT / "docs/parity-motion.json"
EVIDENCE = {
    "independent_motion", "stationary_force", "qualified_motion",
    "mixed_modes", "open_motion",
}


def rows(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def render():
    cases = rows(INVENTORY)
    statuses = rows(STATUS)
    inventory_names = [row["case"] for row in cases]
    status_names = [row["case"] for row in statuses]
    if len(set(inventory_names)) != len(cases):
        raise ValueError("duplicate published case in inventory")
    if len(set(status_names)) != len(statuses):
        raise ValueError("duplicate published case in evidence ledger")
    if set(inventory_names) != set(status_names):
        missing = sorted(set(inventory_names) - set(status_names))
        extra = sorted(set(status_names) - set(inventory_names))
        raise ValueError(f"evidence ledger differs from inventory: missing={missing}, extra={extra}")
    unknown = sorted({row["motion_evidence"] for row in statuses} - EVIDENCE)
    if unknown:
        raise ValueError(f"unknown motion evidence categories: {unknown}")
    gated = sum(row["motion_evidence"] == "independent_motion" for row in statuses)
    payload = {
        "schemaVersion": 1,
        "label": "WEC-Sim case motion",
        "message": f"{gated}/{len(cases)} gated",
        "color": "blue",
    }
    return json.dumps(payload, indent=2) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="fail if the checked-in badge is stale")
    args = parser.parse_args()
    contents = render()
    if args.check:
        if not BADGE.exists() or BADGE.read_text() != contents:
            parser.error("badge is stale; run python tools/render_parity_badge.py")
    else:
        BADGE.parent.mkdir(parents=True, exist_ok=True)
        BADGE.write_text(contents)


if __name__ == "__main__":
    main()

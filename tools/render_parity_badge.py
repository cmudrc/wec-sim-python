"""Render disjoint reference-model and additional-system parity badges."""

import argparse
import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
INVENTORY = ROOT / "tests/reference_cases.csv"
STATUS = ROOT / "tests/reference_case_status.csv"
REFERENCE_FAMILIES = {"RM3", "OSWEC", "Sphere"}
ADDITIONAL_FAMILY = "Other or dynamic"
BADGES = (
    (ROOT / "docs/reference-model-cases.json", "RM3 / OSWEC / Sphere", REFERENCE_FAMILIES),
    (ROOT / "docs/additional-system-cases.json", "Additional systems", {ADDITIONAL_FAMILY}),
)
EVIDENCE = {
    "independent_motion", "stationary_force", "qualified_motion",
    "mixed_modes", "open_motion",
}
GATED = {"independent_motion", "stationary_force"}


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
    status_by_case = {row["case"]: row["motion_evidence"] for row in statuses}
    if {row["family"] for row in cases} != REFERENCE_FAMILIES | {ADDITIONAL_FAMILY}:
        raise ValueError("inventory families differ from badge groups")
    rendered = {}
    for path, label, families in BADGES:
        names = [row["case"] for row in cases if row["family"] in families]
        gated = sum(status_by_case[name] in GATED for name in names)
        payload = {
            "schemaVersion": 1,
            "label": label,
            "message": f"{gated}/{len(names)} gated",
            "color": "blue",
        }
        rendered[path] = json.dumps(payload, indent=2) + "\n"
    return rendered


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="fail if the checked-in badge is stale")
    args = parser.parse_args()
    badges = render()
    if args.check:
        for path, contents in badges.items():
            if not path.exists() or path.read_text() != contents:
                parser.error(f"{path.name} is stale; run python tools/render_parity_badge.py")
    else:
        for path, contents in badges.items():
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(contents)


if __name__ == "__main__":
    main()

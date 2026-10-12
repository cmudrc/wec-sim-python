"""Render the Shields badge from pytest's combined line and branch coverage."""

import argparse
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
COVERAGE = ROOT / "coverage.json"
BADGE = ROOT / "docs/coverage.json"


def render():
    totals = json.loads(COVERAGE.read_text())["totals"]
    percent = int(totals["percent_covered_display"])
    if not 0 <= percent <= 100:
        raise ValueError("coverage percentage is outside 0–100")
    color = ("red" if percent < 50 else "orange" if percent < 70
             else "yellow" if percent < 80 else "green")
    return json.dumps({
        "schemaVersion": 1,
        "label": "Test coverage",
        "message": f"{percent}%",
        "color": color,
    }, indent=2) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true",
                        help="fail if the checked-in badge is stale")
    args = parser.parse_args()
    contents = render()
    if args.check:
        if not BADGE.exists() or BADGE.read_text() != contents:
            parser.error("badge is stale; run python tools/render_coverage_badge.py")
    else:
        BADGE.write_text(contents)


if __name__ == "__main__":
    main()

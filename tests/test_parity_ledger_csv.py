"""Keep published-case evidence intact when rendering the README badges."""

import pytest

from tools.render_parity_badge import STATUS, render, rows


def test_published_case_status_keeps_the_full_wavestar_limitation():
    status = {row["case"]: row for row in rows(STATUS)}
    gap = status["WECCCOMP/WECCCOMP_Nonlinear_Model_Predictive"]["remaining_gap"]
    assert "0.2392% the largest relative error" in gap
    assert "strict badge count is unchanged" in gap
    render()


@pytest.mark.parametrize("record", (
    "WaveStar,qualified_motion,first part,missing quote\n",
    "WaveStar,qualified_motion\n",
))
def test_badge_input_rejects_incomplete_or_overflowing_evidence(tmp_path, record):
    path = tmp_path / "bad_status.csv"
    path.write_text("case,motion_evidence,remaining_gap\n" + record)
    with pytest.raises(ValueError, match="invalid column count"):
        rows(path)

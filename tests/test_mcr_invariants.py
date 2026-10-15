"""Portable multiple-condition input and signed-power behavior."""

import numpy as np
import pytest
from openpyxl import Workbook
from scipy.io import savemat

from wecsim.analysis.mcr import (
    MCRCondition, MCRTrace, mcr_grid, mcr_mat_file, mcr_wave_statistics,
    run_mcr,
)


def test_mcr_condition_order_and_signed_average_power_matrix():
    cases = mcr_grid([1, 2], [6, 8], [100, 200])
    assert [(case.wave_height, case.wave_period, case.pto_damping)
            for case in cases[:4]] == [
                (1, 6, 100), (1, 8, 100), (2, 6, 100), (2, 8, 100),
            ]

    def simulate(case):
        scale = case.wave_height * case.wave_period * case.pto_damping / 100
        return MCRTrace(np.array([0., 1., 2.]),
                        scale * np.array([-5., 2., 4.]))

    result = run_mcr(cases, simulate, averaging_start_time=1)
    assert len(result.traces) == 8
    np.testing.assert_allclose(result.mean_absorbed_power[:4], [18, 24, 36, 48])
    matrix = result.power_matrix(damping=100)
    np.testing.assert_array_equal(matrix.periods, [8, 6])
    np.testing.assert_array_equal(matrix.heights, [1, 2])
    np.testing.assert_array_equal(matrix.absorbed_power, [[24, 48], [18, 36]])
    with pytest.raises(ValueError, match="no MCR conditions"):
        result.power_matrix(damping=999)


def test_mcr_excel_style_selection_and_mat_table_have_same_case_order(tmp_path):
    grid = np.array([[0, 6, 8], [1, 1, 0], [2, 0, 1]])
    selected = mcr_wave_statistics(grid, [100, 200])
    assert selected == (
        MCRCondition(1, 6, 100), MCRCondition(2, 8, 100),
        MCRCondition(1, 6, 200), MCRCondition(2, 8, 200),
    )
    mat = tmp_path / "mcr.mat"
    savemat(mat, {"mcr": {
        "header": np.array(["waves.height", "waves.period",
                            "pto(1).damping", "pto(1).stiffness"], dtype=object),
        "cases": np.array([[1, 6, 100, 0], [2, 8, 100, 0],
                           [1, 6, 200, 0], [2, 8, 200, 0]]),
    }})
    assert mcr_mat_file(mat) == selected
    with pytest.raises(ValueError, match="numbers"):
        mcr_wave_statistics([[0, 6], [1, "=2+2"]], [100])


def test_mcr_excel_workbook_preserves_sea_state_order_and_rejects_formulas(tmp_path):
    workbook = Workbook()
    sheet = workbook.active
    for row in ((None, 6, 8), (1.5, 3, 0), (2.5, 0, 7)):
        sheet.append(row)
    path = tmp_path / "waveStatistic.xlsx"
    workbook.save(path)

    # Positive interior cells select cases; their magnitudes are not weights.
    assert mcr_wave_statistics(path, [1_200_000, 2_400_000]) == (
        MCRCondition(1.5, 6, 1_200_000),
        MCRCondition(2.5, 8, 1_200_000),
        MCRCondition(1.5, 6, 2_400_000),
        MCRCondition(2.5, 8, 2_400_000),
    )

    # A formula has no stable value without Excel recalculation. Reading the
    # formula itself prevents a silently cached result from selecting a case.
    sheet["B2"] = "=1+2"
    workbook.save(path)
    workbook.close()
    with pytest.raises(ValueError, match="numbers, not formulas"):
        mcr_wave_statistics(path, [1_200_000])


def test_mcr_rejects_missing_or_duplicate_power_matrix_cells():
    time = np.array([0., 1.])
    trace = lambda _case: MCRTrace(time, np.array([0., 1.]))
    duplicate = run_mcr([MCRCondition(1, 6, 100)] * 2, trace,
                        averaging_start_time=0)
    with pytest.raises(ValueError, match="duplicate"):
        duplicate.power_matrix(damping=100)
    missing = run_mcr([MCRCondition(1, 6, 100), MCRCondition(2, 8, 100)],
                      trace, averaging_start_time=0)
    with pytest.raises(ValueError, match="incomplete"):
        missing.power_matrix(damping=100)

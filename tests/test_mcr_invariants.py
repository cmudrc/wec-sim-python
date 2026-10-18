"""Portable multiple-condition input and signed-power behavior."""

from pathlib import Path
from shutil import copyfile

import numpy as np
import pytest
from openpyxl import Workbook
from scipy.io import savemat

from wecsim.analysis.mcr import (
    MCRCondition, MCRTrace, mcr_grid, mcr_mat_file, mcr_wave_statistics,
    mcr_spectrum_files, run_mcr, run_rm3_spectrum_mcr,
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


@pytest.mark.parametrize("record", (
    {"cases": np.array([[1, 6, 100, 0]])},
    {"header": np.array(["waves.height", "waves.period",
                         "pto(1).damping", "pto(1).stiffness"], dtype=object)},
))
def test_mcr_mat_table_reports_missing_fields_as_input_error(tmp_path, record):
    path = tmp_path / "incomplete-mcr.mat"
    savemat(path, {"mcr": record})
    with pytest.raises(ValueError, match="MCR MAT file needs header and cases"):
        mcr_mat_file(path)


def test_single_imported_sea_runs_and_reports_its_own_power(tmp_path):
    """MATLAB's one-row cell table still represents one complete MCR case."""
    spectrum = tmp_path / "sea.mat"
    savemat(spectrum, {"spectrumData": np.array([
        [.10, .20, .1], [.125, .35, .4], [.15, .20, .7],
    ])})
    table = tmp_path / "mcr.mat"
    savemat(table, {"mcr": {
        "header": np.array(["waves.spectrumFile", "simu.solver"],
                           dtype=object),
        "cases": np.array([[spectrum.name, "ode4"]], dtype=object),
    }})
    assert mcr_spectrum_files(table) == (spectrum.resolve(),)

    # Use the bundled HDF5 from the repository, independent of the MAT table.
    hydro = Path(__file__).resolve().parents[1] / "examples/data/rm3.h5"
    result = run_rm3_spectrum_mcr(
        hydro, table, dt=.05, end_time=.2, ramp_time=0,
        radiation_memory=.1, averaging_start_time=.1,
    )
    assert result.spectrum_files == (spectrum.resolve(),)
    assert len(result.traces) == len(result.wave_elevation) == 1
    trace = result.traces[0]
    np.testing.assert_array_equal(trace.time, np.arange(5) * .05)
    assert np.ptp(result.wave_elevation[0]) > 0
    np.testing.assert_allclose(
        result.mean_absorbed_power[0],
        np.mean(trace.absorbed_power[2:]), rtol=0, atol=1e-12,
    )

    second_spectrum = tmp_path / "second_sea.mat"
    savemat(second_spectrum, {"spectrumData": np.array([
        [.10, .20, 1.1], [.125, .35, 1.4], [.15, .20, 1.7],
    ])})
    two_seas = tmp_path / "two_seas.mat"
    savemat(two_seas, {"mcr": {
        "header": np.array(["waves.spectrumFile", "simu.solver"],
                           dtype=object),
        "cases": np.array([[spectrum.name, "ode4"],
                           [second_spectrum.name, "ode4"]], dtype=object),
    }})
    combined = run_rm3_spectrum_mcr(
        hydro, two_seas, dt=.05, end_time=.2, ramp_time=0,
        radiation_memory=.1, averaging_start_time=.1,
    )
    assert combined.spectrum_files == (spectrum.resolve(),
                                       second_spectrum.resolve())
    np.testing.assert_array_equal(combined.wave_elevation[0],
                                  result.wave_elevation[0])
    np.testing.assert_array_equal(combined.traces[0].absorbed_power,
                                  trace.absorbed_power)
    assert np.ptp(combined.wave_elevation[0] - combined.wave_elevation[1]) > 0


def test_two_column_mcr_seed_matches_local_matlab_sea(tmp_path):
    """MCR must pass its requested phase realization through to the solver."""
    root = Path(__file__).resolve().parents[1]
    source = root / "tests/fixtures/matlab/imported_spectrum_two_column"
    copyfile(source / "spectrum.mat", tmp_path / "sea.mat")
    table = tmp_path / "mcr.mat"
    savemat(table, {"mcr": {
        "header": np.array(["waves.spectrumFile", "simu.solver"],
                           dtype=object),
        "cases": np.array([["sea.mat", "ode4"]], dtype=object),
    }})
    hydro = root / "examples/data/rm3.h5"
    settings = dict(dt=.1, end_time=2, ramp_time=1,
                    radiation_memory=1, averaging_start_time=1)
    first = run_rm3_spectrum_mcr(
        hydro, table, seed=1, phase_generator="matlab", **settings,
    )
    matlab = np.loadtxt(source / "elevation.csv", delimiter=",")
    np.testing.assert_allclose(first.wave_elevation[0], matlab[:21, 1],
                               rtol=0, atol=1e-12)
    assert np.isfinite(first.traces[0].absorbed_power).all()

    second = run_rm3_spectrum_mcr(
        hydro, table, seed=2, phase_generator="matlab", **settings,
    )
    assert not np.array_equal(second.wave_elevation[0],
                              first.wave_elevation[0])


def test_imported_sea_table_requires_case_rows_and_two_columns(tmp_path):
    """A transposed cell table must not masquerade as one valid sea state."""
    table = tmp_path / "transposed.mat"
    savemat(table, {"mcr": {
        "header": np.array(["waves.spectrumFile", "simu.solver"],
                           dtype=object),
        "cases": np.array([["sea.mat"], ["ode4"]], dtype=object),
    }})
    with pytest.raises(ValueError, match="spectrum and solver columns"):
        mcr_spectrum_files(table)


@pytest.mark.parametrize("header, cases, message", [
    (["waves.spectrumFile", "simu.solver"], [["sea.mat", "ode45"]],
     "fixed-step ode4"),
    (["waves.spectrumFile", "simu.solver"], [["../sea.mat", "ode4"]],
     "local to the MAT table"),
    (["waves.spectrumFile", "simu.solver"],
     [["sea.mat", "ode4"], ["sea.mat", "ode4"]], "duplicated"),
    (["waves.spectrumFile", "simu.step"], [["sea.mat", "ode4"]],
     "unsupported imported-spectrum MCR fields"),
])
def test_imported_sea_table_rejects_cases_that_change_published_inputs(
        tmp_path, header, cases, message):
    """MCR must not silently substitute a solver or leave its input folder."""
    table = tmp_path / "invalid-mcr.mat"
    savemat(table, {"mcr": {
        "header": np.asarray(header, dtype=object),
        "cases": np.asarray(cases, dtype=object),
    }})
    with pytest.raises(ValueError, match=message):
        mcr_spectrum_files(table)


@pytest.mark.parametrize("time, power, averaging_start, message", [
    ([0, 1, 1], [0, 1, 2], 0, "increasing time"),
    ([0, 1, 2], [0, np.nan, 2], 0, "finite absorbed power"),
    ([0, 1, 2], [0, 1], 0, "finite absorbed power"),
    ([0, 1, 2], [0, 1, 2], 3, "beyond an MCR trace"),
])
def test_mcr_rejects_traces_that_would_corrupt_averaged_power(
        time, power, averaging_start, message):
    condition = MCRCondition(1, 6, 100)
    simulate = lambda _: MCRTrace(np.asarray(time), np.asarray(power))
    with pytest.raises(ValueError, match=message):
        run_mcr([condition], simulate, averaging_start_time=averaging_start)


def test_mcr_rejects_nonphysical_settings_and_missing_trace():
    with pytest.raises(ValueError, match="physical wave and PTO settings"):
        MCRCondition(1, 0, 100)
    with pytest.raises(TypeError, match="MCRTrace"):
        run_mcr([MCRCondition(1, 6, 100)], lambda _: None,
                averaging_start_time=0)

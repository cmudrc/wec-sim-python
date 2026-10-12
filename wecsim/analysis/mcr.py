"""Python multiple-condition runs and PTO power matrices.

The case order follows the pinned WEC-Sim ``wecSimMCR.m`` script: sea states
vary fastest, then PTO damping, then PTO stiffness. The published RM3 MAT
file and Excel workbook are supported input formats. Simulation itself runs
through the existing Python dynamics solver.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable

import numpy as np
from scipy.io import loadmat

from wecsim.dynamics.rm3Regular import solve_rm3_regular
from wecsim.environment.irregularWave import imported_spectrum_components, synthesize_irregular_response


@dataclass(frozen=True)
class MCRCondition:
    wave_height: float
    wave_period: float
    pto_damping: float
    pto_stiffness: float = 0.0

    def __post_init__(self):
        values = (self.wave_height, self.wave_period,
                  self.pto_damping, self.pto_stiffness)
        if (not np.isfinite(values).all() or self.wave_height < 0
                or self.wave_period <= 0 or self.pto_damping < 0
                or self.pto_stiffness < 0):
            raise ValueError("MCR condition needs finite, physical wave and PTO settings")


@dataclass(frozen=True)
class MCRTrace:
    time: np.ndarray
    absorbed_power: np.ndarray
    response: object | None = None


@dataclass(frozen=True)
class MCRPowerMatrix:
    """Power rows are periods descending; columns are heights ascending."""

    periods: np.ndarray
    heights: np.ndarray
    absorbed_power: np.ndarray


@dataclass(frozen=True)
class MCRResult:
    conditions: tuple[MCRCondition, ...]
    mean_absorbed_power: np.ndarray
    traces: tuple[MCRTrace, ...]

    def power_matrix(self, *, damping: float,
                     stiffness: float = 0.0) -> MCRPowerMatrix:
        """Arrange a complete sea-state grid like the published RM3 figures."""
        selected = [
            (condition, power)
            for condition, power in zip(self.conditions, self.mean_absorbed_power)
            if condition.pto_damping == damping
            and condition.pto_stiffness == stiffness
        ]
        if not selected:
            raise ValueError("no MCR conditions match these PTO settings")
        heights = np.array(sorted({case.wave_height for case, _ in selected}))
        periods = np.array(sorted({case.wave_period for case, _ in selected},
                                  reverse=True))
        matrix = np.full((len(periods), len(heights)), np.nan)
        for condition, power in selected:
            row = int(np.flatnonzero(periods == condition.wave_period)[0])
            column = int(np.flatnonzero(heights == condition.wave_height)[0])
            if np.isfinite(matrix[row, column]):
                raise ValueError("duplicate MCR sea state at these PTO settings")
            matrix[row, column] = power
        if not np.isfinite(matrix).all():
            raise ValueError("MCR sea-state grid is incomplete for these PTO settings")
        return MCRPowerMatrix(periods, heights, matrix)


@dataclass(frozen=True)
class MCRSeaStateResult:
    """Published spectrum-import MCR outputs in source case order."""

    spectrum_files: tuple[Path, ...]
    mean_absorbed_power: np.ndarray
    traces: tuple[MCRTrace, ...]
    wave_elevation: tuple[np.ndarray, ...]


def mcr_spectrum_files(path: str | Path) -> tuple[Path, ...]:
    """Read the published Option 3 sea-state MAT table and resolve its files."""
    mat_file = Path(path).expanduser().resolve(strict=True)
    data = loadmat(mat_file, squeeze_me=True, struct_as_record=False)
    if "mcr" not in data:
        raise ValueError("MAT file has no mcr structure")
    mcr = data["mcr"]
    if (not hasattr(mcr, "header") or not hasattr(mcr, "cases")
            or tuple(str(value) for value in np.ravel(mcr.header))
            != ("waves.spectrumFile", "simu.solver")):
        raise ValueError("unsupported imported-spectrum MCR fields")
    cases = np.asarray(mcr.cases, dtype=object)
    if cases.ndim != 2 or cases.shape[1] != 2 or len(cases) == 0:
        raise ValueError("imported-spectrum MCR needs spectrum and solver columns")
    if any(str(solver) != "ode4" for solver in cases[:, 1]):
        raise ValueError("only the published fixed-step ode4 solver is supported")
    filenames = tuple(str(value) for value in cases[:, 0])
    if any(not name or Path(name).name != name for name in filenames):
        raise ValueError("spectrum file names must be local to the MAT table")
    if len(set(filenames)) != len(filenames):
        raise ValueError("spectrum files must not be duplicated")
    return tuple((mat_file.parent / name).resolve(strict=True)
                 for name in filenames)


def _values(values, name, *, positive=False, nonnegative=False):
    array = np.asarray(tuple(values), dtype=float)
    if (array.ndim != 1 or not len(array) or not np.isfinite(array).all()
            or (positive and np.any(array <= 0))
            or (nonnegative and np.any(array < 0))):
        raise ValueError(f"{name} must contain finite supported values")
    if len(np.unique(array)) != len(array):
        raise ValueError(f"{name} must not contain duplicates")
    return array


def _expand(sea_states, damping_values, stiffness_values):
    damping = _values(damping_values, "PTO damping", nonnegative=True)
    stiffness = _values(stiffness_values, "PTO stiffness", nonnegative=True)
    if not sea_states:
        raise ValueError("MCR needs at least one sea state")
    return tuple(
        MCRCondition(height, period, coefficient, spring)
        for spring in stiffness
        for coefficient in damping
        for height, period in sea_states
    )


def mcr_grid(
    heights: Iterable[float], periods: Iterable[float],
    damping_values: Iterable[float],
    stiffness_values: Iterable[float] = (0.0,),
) -> tuple[MCRCondition, ...]:
    """Expand MATLAB MCR Option 1 arrays in source case order."""
    height = _values(heights, "wave heights", nonnegative=True)
    period = _values(periods, "wave periods", positive=True)
    return _expand([(h, p) for h in height for p in period],
                   damping_values, stiffness_values)


def mcr_wave_statistics(
    matrix_or_file: np.ndarray | str | Path,
    damping_values: Iterable[float],
    stiffness_values: Iterable[float] = (0.0,),
) -> tuple[MCRCondition, ...]:
    """Expand positive sea-state cells in the Option 2 Excel grid.

    The first row contains periods, the first column heights. A positive
    interior value selects that sea state, as in ``wecSimMCR.m``. The value is
    not a weight in the published postprocessing.
    Workbook paths require the ``excel`` install extra; an in-memory grid does not.
    """
    if isinstance(matrix_or_file, (str, Path)):
        from openpyxl import load_workbook

        workbook = load_workbook(matrix_or_file, read_only=True,
                                 data_only=False)
        try:
            rows = list(workbook.active.values)
        finally:
            workbook.close()
    else:
        rows = np.asarray(matrix_or_file, dtype=object).tolist()
    if not isinstance(rows, list) or len(rows) < 2:
        raise ValueError("wave-statistics grid needs a header and sea-state rows")
    width = max(map(len, rows))
    if width < 2:
        raise ValueError("wave-statistics grid needs period columns")
    grid = np.zeros((len(rows), width))
    for row_index, row in enumerate(rows):
        for column_index, value in enumerate(row):
            if value is None:
                continue
            if isinstance(value, (bool, str)):
                raise ValueError("wave-statistics cells must contain numbers, not formulas or text")
            try:
                grid[row_index, column_index] = float(value)
            except (TypeError, ValueError) as exc:
                raise ValueError("wave-statistics cells must be numeric") from exc
    if not np.isfinite(grid).all():
        raise ValueError("wave-statistics cells must be finite")
    periods = grid[0, 1:]
    heights = grid[1:, 0]
    if (np.any(periods <= 0) or np.any(heights < 0)
            or len(np.unique(periods)) != len(periods)
            or len(np.unique(heights)) != len(heights)):
        raise ValueError("wave-statistics heights and periods must be valid and unique")
    sea_states = [
        (heights[i], periods[j])
        for i in range(len(heights))
        for j in range(len(periods))
        if grid[i + 1, j + 1] > 0
    ]
    return _expand(sea_states, damping_values, stiffness_values)


def mcr_mat_file(path: str | Path) -> tuple[MCRCondition, ...]:
    """Read the published Option 3 MAT-file condition table."""
    data = loadmat(path, squeeze_me=True, struct_as_record=False)
    if "mcr" not in data:
        raise ValueError("MAT file has no mcr structure")
    mcr = data["mcr"]
    header = tuple(str(value) for value in np.ravel(mcr.header))
    if header != ("waves.height", "waves.period", "pto(1).damping",
                  "pto(1).stiffness"):
        raise ValueError("unsupported MCR MAT-file fields")
    cases = np.atleast_2d(np.asarray(mcr.cases, dtype=float))
    if cases.shape[1] != 4 or not np.isfinite(cases).all():
        raise ValueError("MCR MAT file needs four finite columns")
    return tuple(MCRCondition(*row) for row in cases)


def run_mcr(
    conditions: Iterable[MCRCondition],
    simulate: Callable[[MCRCondition], MCRTrace],
    *,
    averaging_start_time: float,
) -> MCRResult:
    """Run each configured WEC and average signed absorbed PTO power.

    A spring can return energy during part of a cycle, so instantaneous
    absorbed power may be negative. The averaging window includes its first
    sample, matching MATLAB's ``mean(power(2000:end))`` at 0.1 s steps.
    """
    cases = tuple(conditions)
    if not cases or not all(isinstance(case, MCRCondition) for case in cases):
        raise ValueError("run_mcr needs MCRCondition values")
    if not np.isfinite(averaging_start_time) or averaging_start_time < 0:
        raise ValueError("averaging start time must be finite and nonnegative")
    traces = []
    averages = []
    for condition in cases:
        trace = simulate(condition)
        if not isinstance(trace, MCRTrace):
            raise TypeError("simulate must return MCRTrace")
        time = np.asarray(trace.time, dtype=float)
        power = np.asarray(trace.absorbed_power, dtype=float)
        if (time.ndim != 1 or len(time) < 2 or power.shape != time.shape
                or not np.isfinite(time).all() or not np.isfinite(power).all()
                or time[0] < 0 or np.any(np.diff(time) <= 0)):
            raise ValueError("MCR trace needs increasing time and finite absorbed power")
        start = int(np.searchsorted(time, averaging_start_time - 1e-10))
        if start >= len(time):
            raise ValueError("averaging start time is beyond an MCR trace")
        traces.append(MCRTrace(time.copy(), power.copy(), trace.response))
        averages.append(float(np.mean(power[start:])))
    return MCRResult(cases, np.asarray(averages), tuple(traces))


def run_rm3_mcr(
    h5_file: str | Path, conditions: Iterable[MCRCondition], *,
    dt: float = 0.1, end_time: float = 400.0,
    ramp_time: float = 100.0, radiation_memory: float = 60.0,
    averaging_start_time: float = 199.9,
    added_mass_scheme: str = "implicit",
) -> MCRResult:
    """Run the published RM3 floating-joint conditions and PTO power matrix."""
    def simulate(condition):
        response = solve_rm3_regular(
            h5_file, wave_height=condition.wave_height,
            wave_period=condition.wave_period,
            pto_damping=condition.pto_damping,
            pto_stiffness=condition.pto_stiffness,
            radiation_memory=radiation_memory,
            added_mass_scheme=added_mass_scheme,
            dt=dt, end_time=end_time, ramp_time=ramp_time,
        )
        return MCRTrace(response.time, response.pto_mechanical_power, response)

    return run_mcr(conditions, simulate,
                   averaging_start_time=averaging_start_time)


def run_rm3_spectrum_mcr(
    h5_file: str | Path, mat_file: str | Path, *,
    dt: float = 0.1, end_time: float = 400.0,
    ramp_time: float = 100.0, radiation_memory: float = 60.0,
    pto_damping: float = 1_200_000.0,
    averaging_start_time: float = 199.9,
    added_mass_scheme: str = "implicit",
) -> MCRSeaStateResult:
    """Run the published three imported-spectrum RM3 MCR sea states.

    Spectrum frequencies, densities, and phases come from each source MAT
    file. The three cases share the published RM3 floating joint and PTO.
    """
    if (not np.isfinite(averaging_start_time)
            or averaging_start_time < 0
            or averaging_start_time > end_time):
        raise ValueError("averaging start time must be within the run")
    files = mcr_spectrum_files(mat_file)
    traces = []
    means = []
    elevations = []
    for spectrum_file in files:
        components = imported_spectrum_components(h5_file, spectrum_file)
        wave = tuple(
            synthesize_irregular_response(
                h5_file, components, dt=dt, end_time=end_time,
                ramp_time=ramp_time, body_number=body,
            )
            for body in (1, 2)
        )
        forcing = np.stack([response.excitation_force for response in wave], axis=1)
        response = solve_rm3_regular(
            h5_file, wave_height=0, pto_damping=pto_damping,
            radiation_memory=radiation_memory,
            added_mass_scheme=added_mass_scheme,
            excitation_force=forcing,
            dt=dt, end_time=end_time, ramp_time=ramp_time,
        )
        start = int(np.searchsorted(response.time,
                                    averaging_start_time - 1e-10))
        power = response.pto_mechanical_power.copy()
        traces.append(MCRTrace(response.time.copy(), power, response))
        means.append(float(np.mean(power[start:])))
        elevations.append(wave[0].elevation.copy())
    return MCRSeaStateResult(files, np.asarray(means),
                             tuple(traces), tuple(elevations))

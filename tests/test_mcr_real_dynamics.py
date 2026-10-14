"""Short real RM3 MCR runs that need no MATLAB installation or source files."""

from pathlib import Path

import numpy as np
from scipy.io import savemat

from wecsim import mcr_grid, run_rm3_mcr, run_rm3_spectrum_mcr


HYDRO = Path(__file__).resolve().parents[1] / "examples/data/rm3.h5"


def test_regular_mcr_averages_real_spring_damper_power_in_case_order():
    conditions = mcr_grid([1.5, 2.5], [6, 8], [1_200_000, 2_400_000],
                          stiffness_values=[50_000])
    result = run_rm3_mcr(
        HYDRO, conditions, dt=.1, end_time=1, ramp_time=0,
        radiation_memory=.2, averaging_start_time=.5,
    )
    assert result.conditions == conditions
    for condition, trace, mean in zip(
            conditions, result.traces, result.mean_absorbed_power):
        response = trace.response
        speed, stroke = response.pto_velocity, response.pto_stroke
        np.testing.assert_allclose(
            response.pto_force,
            -condition.pto_damping * speed - condition.pto_stiffness * stroke,
            rtol=0, atol=1e-7,
        )
        np.testing.assert_allclose(trace.absorbed_power,
                                   -response.pto_force * speed,
                                   rtol=0, atol=1e-7)
        np.testing.assert_allclose(response.pto_dissipated_power,
                                   condition.pto_damping * speed**2,
                                   rtol=0, atol=1e-7)
        # The spring stores energy, so mechanical and dissipated power differ.
        assert np.max(np.abs(trace.absorbed_power
                             - response.pto_dissipated_power)) > 10
        np.testing.assert_allclose(mean,
                                   trace.absorbed_power[trace.time >= .5].mean(),
                                   rtol=0, atol=1e-10)

    for damping, offset in ((1_200_000, 0), (2_400_000, 4)):
        matrix = result.power_matrix(damping=damping, stiffness=50_000)
        np.testing.assert_array_equal(matrix.periods, [8, 6])
        np.testing.assert_array_equal(matrix.heights, [1.5, 2.5])
        np.testing.assert_array_equal(
            matrix.absorbed_power,
            [[result.mean_absorbed_power[offset + 1],
              result.mean_absorbed_power[offset + 3]],
             [result.mean_absorbed_power[offset],
              result.mean_absorbed_power[offset + 2]]],
        )
        assert np.all(matrix.absorbed_power[:, 1] > matrix.absorbed_power[:, 0])
        assert np.all(matrix.absorbed_power > 0)
        assert np.all(matrix.absorbed_power[0] > matrix.absorbed_power[1])
    assert result.mean_absorbed_power[1] != result.mean_absorbed_power[5]


def test_imported_spectrum_mcr_uses_each_sea_and_averages_actual_pto_power(tmp_path):
    frequency = np.array([.1, .15, .2, .25])
    density = np.array([.2, .4, .3, .1])
    phase = np.array([0, .2, .4, .6])
    for index, scale in enumerate((1., 2.), start=1):
        savemat(tmp_path / f"sea{index}.mat", {
            "spectrumData": np.column_stack((frequency, scale * density, phase)),
        })
    table = tmp_path / "mcr.mat"
    savemat(table, {"mcr": {
        "header": np.array(["waves.spectrumFile", "simu.solver"], dtype=object),
        "cases": np.array([["sea1.mat", "ode4"],
                           ["sea2.mat", "ode4"]], dtype=object),
    }})
    result = run_rm3_spectrum_mcr(
        HYDRO, table, dt=.1, end_time=1, ramp_time=0,
        radiation_memory=.2, averaging_start_time=.5,
    )
    assert tuple(path.name for path in result.spectrum_files) == (
        "sea1.mat", "sea2.mat",
    )
    # Doubling spectral density raises incident wave amplitude by sqrt(2).
    np.testing.assert_allclose(result.wave_elevation[1],
                               np.sqrt(2) * result.wave_elevation[0],
                               rtol=0, atol=1e-12)
    for trace, mean in zip(result.traces, result.mean_absorbed_power):
        response = trace.response
        np.testing.assert_allclose(trace.absorbed_power,
                                   1_200_000 * response.pto_velocity**2,
                                   rtol=0, atol=1e-7)
        np.testing.assert_allclose(mean,
                                   trace.absorbed_power[trace.time >= .5].mean(),
                                   rtol=0, atol=1e-10)
    assert result.mean_absorbed_power[1] > result.mean_absorbed_power[0] > 0

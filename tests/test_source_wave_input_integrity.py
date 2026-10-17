"""Guard the MATLAB MAT sea inputs used by paired motion comparisons."""

from pathlib import Path

import numpy as np
import pytest
from scipy.io import savemat

from wecsim.environment.irregularWave import (
    imported_full_directional_components,
    imported_spectrum_components,
    synthesize_irregular_response,
)


ROOT = Path(__file__).resolve().parents[1]
OSWEC = ROOT / "examples/data/oswec.h5"
SPHERE = ROOT / "examples/data/sphere.h5"


def _full_directional_data():
    width = np.deg2rad(30)
    return {
        "frequencies": np.array([.1, .15, .2]),
        "spectrum": np.array([.2, .4, .2]),
        "directions": np.array([0., 30., 60.]),
        "spread": np.full((3, 3), 1 / (3 * width)),
    }


@pytest.mark.parametrize("damage, message", (
    ("missing_spread", "needs frequencies"),
    ("negative_density", "invalid frequency"),
    ("repeated_heading", "invalid frequency"),
    ("one_bem_bin", "two BEM-range frequencies"),
))
def test_full_directional_import_rejects_corrupt_source_arrays(tmp_path,
                                                                damage, message):
    data = _full_directional_data()
    if damage == "missing_spread":
        del data["spread"]
    elif damage == "negative_density":
        data["spectrum"][1] = -.4
    elif damage == "repeated_heading":
        data["directions"][1] = 0
    else:
        data["frequencies"] = np.array([.0001, .0002, .1])
    path = tmp_path / "full_directional.mat"
    savemat(path, data)
    with pytest.raises(ValueError, match=message):
        imported_full_directional_components(OSWEC, path, phase=np.zeros((3, 3)))


def test_full_directional_import_requires_phase_grid_of_retained_bins(tmp_path):
    path = tmp_path / "full_directional.mat"
    savemat(path, _full_directional_data())
    with pytest.raises(ValueError, match="either phase or seed"):
        imported_full_directional_components(
            OSWEC, path, phase=np.zeros((3, 3)), seed=1,
        )
    with pytest.raises(ValueError, match="selected frequency-by-heading grid"):
        imported_full_directional_components(OSWEC, path, phase=np.zeros((2, 3)))


def test_imported_spectrum_keeps_saved_phases_after_bem_filtering(tmp_path):
    # Out-of-band bins must not shift the saved phases of the three retained
    # components: those phases identify a particular MATLAB sea realization.
    values = np.array([
        [.0001, 9., .7], [.1, .2, .1], [.15, .4, 1.2],
        [.2, .2, 2.1], [5., 9., .8],
    ])
    path = tmp_path / "spectrum.mat"
    savemat(path, {"spectrumData": values})
    components = imported_spectrum_components(SPHERE, path)
    np.testing.assert_array_equal(components.omega, 2 * np.pi * values[1:4, 0])
    np.testing.assert_array_equal(components.phase[:, 0], values[1:4, 2])
    response = synthesize_irregular_response(
        SPHERE, components, dt=.1, end_time=.2, ramp_time=0,
    )
    expected = np.sum(np.sqrt(components.spectral_amplitude
                              * components.d_omega)
                      * np.cos(components.phase[:, 0]))
    np.testing.assert_allclose(response.elevation[0], expected,
                               rtol=0, atol=1e-12)
    assert np.isfinite(response.excitation_force).all()


@pytest.mark.parametrize("damage, message", (
    ("missing", "needs spectrumData"),
    ("one_column", "two or three finite columns"),
    ("four_columns", "two or three finite columns"),
    ("nonfinite_phase", "two or three finite columns"),
    ("negative_density", "nonnegative density"),
    ("repeated_frequency", "increasing frequencies"),
    ("one_bem_bin", "increasing frequencies"),
))
def test_imported_spectrum_rejects_corrupt_matlab_bins(tmp_path, damage, message):
    values = np.array([[.1, .2, .1], [.15, .4, 1.2], [.2, .2, 2.1]])
    if damage == "one_column":
        values = values[:, :1]
    elif damage == "four_columns":
        values = np.column_stack((values, np.zeros(3)))
    elif damage == "nonfinite_phase":
        values[1, 2] = np.nan
    elif damage == "negative_density":
        values[1, 1] = -.4
    elif damage == "repeated_frequency":
        values[1, 0] = values[0, 0]
    elif damage == "one_bem_bin":
        values[:2, 0] = [.0001, .0002]
    path = tmp_path / "spectrum.mat"
    savemat(path, {} if damage == "missing" else {"spectrumData": values})
    with pytest.raises(ValueError, match=message):
        imported_spectrum_components(SPHERE, path)

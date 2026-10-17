"""Exercise measured-elevation timing and force ramps through the public API."""

from pathlib import Path

import numpy as np
import pytest
from scipy.io import savemat

from wecsim import ImportedElevationWave, WEC


SPHERE = Path(__file__).resolve().parents[1] / "examples/data/sphere.h5"


def _run(path, *, second_force_ramp=False, variable="etaData", file_name=None):
    wec = WEC("sampled Sphere sea")
    sphere = wec.body("sphere", SPHERE)
    wec.coordinate("heave", sphere.move("heave"))
    return wec.run(
        ImportedElevationWave(path.name if file_name is None else file_name,
                              variable=variable,
                              reapply_force_ramp=second_force_ramp),
        dt=.01, end_time=.2, ramp_time=.1, radiation_memory=.2,
        base_dir=path.parent,
    )


def test_nonuniform_elevation_samples_keep_their_timing_and_force_ramp(tmp_path):
    samples = np.array([
        [0, 0], [.07, 1.3], [.13, -.4], [.23, .6],
    ])
    path = tmp_path / "measured.mat"
    savemat(path, {"etaData": samples})

    once = _run(path)
    twice = _run(path, second_force_ramp=True)
    time = once.time
    ramp = np.ones_like(time)
    early = time < .1
    ramp[early] = (1 - np.cos(np.pi * time[early] / .1)) / 2

    np.testing.assert_array_equal(once.time, twice.time)
    np.testing.assert_allclose(
        once.wave_elevation, np.interp(time, samples[:, 0], samples[:, 1]) * ramp,
        rtol=0, atol=1e-14,
    )
    np.testing.assert_array_equal(once.wave_elevation, twice.wave_elevation)
    first_force = dict(once.raw.extra_outputs)["body1_excitation_force"]
    second_force = dict(twice.raw.extra_outputs)["body1_excitation_force"]
    assert np.max(np.abs(first_force[early])) > 1
    np.testing.assert_allclose(second_force, first_force * ramp[:, None],
                               rtol=0, atol=1e-7)
    assert once.raw.auxiliary_files == twice.raw.auxiliary_files == (path.resolve(),)


@pytest.mark.parametrize(
    ("samples", "message"),
    [
        ([(0, 0), (.1, 1), (.1, 2), (.3, 0)], "finite, increasing"),
        ([(0, 0), (.1, 1), (.2, float("nan")), (.3, 0)], "finite, increasing"),
        ([(0, 0), (.05, 1), (.1, 2), (.15, 0)], "cover the full simulation time"),
    ],
)
def test_imported_elevation_rejects_ambiguous_or_short_records(
    tmp_path, samples, message,
):
    path = tmp_path / "measured.mat"
    savemat(path, {"etaData": np.asarray(samples)})
    with pytest.raises(ValueError, match=message):
        _run(path)


def test_imported_elevation_rejects_missing_or_nonphysical_channels(tmp_path):
    path = tmp_path / "measured.mat"
    samples = np.array([[0, 0], [.1, 1], [.2, 0]])
    savemat(path, {"etaData": samples})
    with pytest.raises(ValueError, match="nonempty MAT variable name"):
        _run(path, variable="")
    with pytest.raises(ValueError, match="is absent"):
        _run(path, variable="wrong_channel")
    with pytest.raises(ValueError, match="must be a boolean"):
        _run(path, second_force_ramp="yes")

    complex_samples = samples.astype(complex)
    complex_samples[1, 1] += 1j
    savemat(path, {"etaData": complex_samples})
    with pytest.raises(ValueError, match="real numeric samples"):
        _run(path)


def test_imported_elevation_rejects_paths_that_are_not_mat_files(tmp_path):
    path = tmp_path / "measured.mat"
    savemat(path, {"etaData": np.array([[0, 0], [.2, 1]])})
    with pytest.raises(ValueError, match="MAT file path"):
        _run(path, file_name="")
    with pytest.raises(ValueError, match="is not a file"):
        _run(tmp_path)

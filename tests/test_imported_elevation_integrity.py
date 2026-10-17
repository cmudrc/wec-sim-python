"""Exercise measured-elevation timing and force ramps through the public API."""

from pathlib import Path

import numpy as np
import pytest
from scipy.io import savemat

from wecsim import ImportedElevationWave, WEC


SPHERE = Path(__file__).resolve().parents[1] / "examples/data/sphere.h5"


def _run(path, *, second_force_ramp=False):
    wec = WEC("sampled Sphere sea")
    sphere = wec.body("sphere", SPHERE)
    wec.coordinate("heave", sphere.move("heave"))
    return wec.run(
        ImportedElevationWave(path.name,
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

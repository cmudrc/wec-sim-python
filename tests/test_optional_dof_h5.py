"""A conventional six-DOF HDF5 file need not repeat its DOF numbering."""

from pathlib import Path
from shutil import copyfile

import h5py
import numpy as np

from wecsim import RegularCICWave, WEC
from wecsim.hydrodynamics.bodyClass import BodyClass


HYDRO = Path(__file__).parent / "fixtures/matlab/sphere_mean_drift/sphere.h5"


def _run(hydro):
    wec = WEC("Sphere optional DOF metadata")
    sphere = wec.body("sphere", hydro, inertia=(837.75804096,) * 3,
                      mean_drift="control_surface")
    for axis in ("surge", "heave", "pitch"):
        wec.coordinate(axis, sphere.move(axis))
    return wec.run(RegularCICWave(0.1, 2), dt=0.01, end_time=2,
                   ramp_time=0, radiation_memory=1).raw


def test_sphere_without_optional_dof_metadata_has_same_short_motion(tmp_path):
    reduced = tmp_path / "sphere-no-dof-metadata.h5"
    copyfile(HYDRO, reduced)
    with h5py.File(reduced, "r+") as file:
        properties = file["body1/properties"]
        for name in ("dof", "dof_start", "dof_end"):
            del properties[name]

    body = BodyClass(str(reduced))
    body.bodyNumber = 1
    body.readH5file()
    np.testing.assert_array_equal(body.dof, [6])
    np.testing.assert_array_equal(body.dof_start, [1])
    np.testing.assert_array_equal(body.dof_end, [6])

    reference = _run(HYDRO)
    actual = _run(reduced)
    assert np.max(np.abs(reference.body_position[:, 0, 2])) > 0
    for name in ("time", "wave_elevation", "body_position", "body_velocity"):
        np.testing.assert_allclose(getattr(actual, name), getattr(reference, name),
                                   rtol=0, atol=1e-12)
    reference_outputs = dict(reference.extra_outputs)
    actual_outputs = dict(actual.extra_outputs)
    assert actual_outputs.keys() == reference_outputs.keys()
    for name in reference_outputs:
        np.testing.assert_allclose(actual_outputs[name], reference_outputs[name],
                                   rtol=0, atol=1e-9)


def test_rm3_optional_dof_defaults_keep_two_body_numbering(tmp_path):
    source = Path(__file__).parents[1] / "examples/data/rm3.h5"
    reduced = tmp_path / "rm3-no-dof-metadata.h5"
    copyfile(source, reduced)
    with h5py.File(reduced, "r+") as file:
        for body_number in (1, 2):
            properties = file[f"body{body_number}/properties"]
            for name in ("dof", "dof_start", "dof_end"):
                del properties[name]

    for body_number, start, end in ((1, 1, 6), (2, 7, 12)):
        original = BodyClass(str(source))
        original.bodyNumber = body_number
        original.readH5file()
        body = BodyClass(str(reduced))
        body.bodyNumber = body_number
        body.readH5file()
        for name, expected in (("dof", 6), ("dof_start", start),
                               ("dof_end", end)):
            np.testing.assert_array_equal(getattr(body, name), [expected])
            np.testing.assert_array_equal(getattr(body, name),
                                          getattr(original, name))

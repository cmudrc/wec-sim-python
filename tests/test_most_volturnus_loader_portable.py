"""Portable source-file and force-law check for the MOST VolturnUS platform."""

import h5py
import numpy as np

from wecsim import MostPlatformHydrodynamics


def test_volturnus_files_produce_expected_restoring_and_radiation(tmp_path):
    """Check the source HDF5 orientation and SI scaling through live forces."""
    rho, gravity = 1025.0, 9.80665
    volume, mass = 10.0, 9000.0
    cg = np.array([0.0, 0.0, -2.0])
    cb = np.array([0.2, 0.0, -1.5])
    inertia = np.diag([1e6, 2e6, 3e6])
    hydro_file = tmp_path / "volturnus.h5"
    mass_file = tmp_path / "mass_properties.h5"

    with h5py.File(hydro_file, "w") as source:
        properties = source.create_group("body1/properties")
        properties.create_dataset("disp_vol", data=[[volume]])
        properties.create_dataset("cg", data=cg[:, None])
        properties.create_dataset("cb", data=cb[:, None])
        coefficients = source.create_group("body1/hydro_coeffs")
        coefficients.create_dataset(
            "linear_restoring_stiffness",
            data=np.diag([0.0, 0.0, 2.0, 0.0, 3.0, 0.0]),
        )
        coefficients.create_dataset(
            "added_mass/inf_freq", data=np.diag([2., 3., 4., 5., 6., 7.]),
        )
        irf = coefficients.create_group("radiation_damping/impulse_response_fun")
        irf.create_dataset("t", data=[0.0, 0.1, 0.2])
        # The source stores lag last; the Python force law needs lag first.
        irf.create_dataset("K", data=np.repeat((0.4 * np.eye(6))[:, :, None], 3,
                                               axis=2))
    with h5py.File(mass_file, "w") as source:
        properties = source.create_group("Platform/VolturnUS")
        properties.create_dataset("mass", data=[[mass]])
        properties.create_dataset("I_COG", data=inertia)

    platform = MostPlatformHydrodynamics.from_volturnus(hydro_file, mass_file)
    assert platform.mass == mass
    np.testing.assert_array_equal(platform.inertia, inertia)
    np.testing.assert_array_equal(platform.equilibrium_pose, np.r_[cg, 0., 0., 0.])
    np.testing.assert_array_equal(platform.added_mass,
                                  rho * np.diag([2., 3., 4., 5., 6., 7.]))

    # At the file's reference pose, buoyancy opposes weight and the offset
    # center of buoyancy creates a pitch moment about the center of gravity.
    buoyancy = rho * gravity * volume
    expected = np.array([0., 0., mass * gravity - buoyancy,
                         0., cb[0] * buoyancy, 0.])
    np.testing.assert_allclose(
        platform.restoring_force(platform.equilibrium_pose), expected,
        rtol=0, atol=1e-10,
    )
    displaced = platform.equilibrium_pose + [0., 0., .3, 0., .02, 0.]
    expected += [0., 0., .3 * 2 * rho * gravity,
                 0., .02 * 3 * rho * gravity, 0.]
    np.testing.assert_allclose(platform.restoring_force(displaced), expected,
                               rtol=0, atol=1e-10)

    # A constant diagonal kernel has an exact trapezoidal convolution on
    # this grid. The first sample has no history; each later lag adds load.
    speed = np.tile([0., 0., .1, 0., 0., 0.], (3, 1))
    radiation = platform.radiation_force(speed, dt=.1, memory_time=.2)
    np.testing.assert_allclose(radiation[:, 2], [0., 4.1, 8.2],
                               rtol=0, atol=1e-12)
    np.testing.assert_array_equal(radiation[:, [0, 1, 3, 4, 5]], 0)

"""Pair the legacy nonlinear body mesh and equilibrium mass with MATLAB data."""

from pathlib import Path

import numpy as np
from scipy.io import loadmat

from wecsim.hydrodynamics.bodyClass import BodyClass


FIXTURE = (Path(__file__).resolve().parent
           / "test_objects/test_bodyclass/testData")
MESH = FIXTURE / "body_10_test"


def _matlab_geometry(field):
    data = loadmat(MESH / f"body10_{field}.mat")
    return next(value for name, value in data.items()
                if not name.startswith("__"))


def test_ellipsoid_mesh_and_nonlinear_equilibrium_mass_match_matlab():
    body = BodyClass(str(FIXTURE / "hydroData/ellipsoid.h5"))
    body.bodyNumber = 1
    body.readH5file()
    body.bodyGeo(MESH / "elipsoid.stl")

    # These facet arrays were exported by MATLAB WEC-Sim for the same STL.
    # Facet orientation matters because submerged volume is a surface integral.
    for field in ("face", "vertex", "center", "area"):
        np.testing.assert_array_equal(body.bodyGeometry[field],
                                      _matlab_geometry(field))
    np.testing.assert_allclose(body.bodyGeometry["norm"],
                               _matlab_geometry("norm"), rtol=0, atol=3e-8)

    density = 1025
    center = _matlab_geometry("center")
    area = _matlab_geometry("area")[:, 0]
    normal_z = _matlab_geometry("norm")[:, 2]
    cg_z = float(body.hydroData["properties"]["cg"][0, 2])
    # MATLAB's nonlinear equilibrium mass integrates only submerged facets.
    submerged_z = np.minimum(center[:, 2] + cg_z, 0)
    source_mass = density * np.sum(submerged_z * area * normal_z)
    body.mass = "equilibrium"
    body.setMassMatrix(density, nlHydro=1)
    np.testing.assert_allclose(body.mass, source_mass, rtol=0, atol=1e-3)

    # The linear path instead uses the HDF5 displaced volume at equilibrium.
    body.mass = "equilibrium"
    body.setMassMatrix(density, nlHydro=0)
    source_volume = float(body.hydroData["properties"]["disp_vol"][0, 0])
    np.testing.assert_allclose(body.mass, density * source_volume,
                               rtol=0, atol=1e-10)

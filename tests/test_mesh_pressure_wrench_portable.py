"""Analytic force and moment checks for moving nonlinear-pressure facets."""

import numpy as np

from wecsim.hydrodynamics.nonlinearHydro import (
    MeshPressures, mesh_pressure_wrenches,
)


def test_pressure_wrenches_rotate_moved_facets_but_not_mean_facets():
    # This triangle has area 4.5 m², centroid (1, 1, -2) m, and +Z normal.
    vertices = np.array([[0., 0., -2.], [3., 0., -2.], [0., 3., -2.]])
    faces = np.array([[0, 1, 2]])
    poses = np.array([
        [0., 0., 0., 0., 0., 0.],
        [0., 0., 0., 0., np.pi / 2, 0.],
    ])
    pressures = MeshPressures(
        hydrostatic=np.full((2, 1), 2.),
        nonlinear_wave=np.full((2, 1), 3.),
        linear_wave=np.full((2, 1), 5.),
    )

    wrench = mesh_pressure_wrenches(vertices, faces, poses, pressures)

    # Force is -pressure × outward area vector; moment is about the CG.
    np.testing.assert_allclose(wrench.hydrostatic,
                               [[0, 0, -9, -9, 9, 0],
                                [-9, 0, 0, 0, 9, 9]], rtol=0, atol=1e-14)
    np.testing.assert_allclose(wrench.nonlinear_wave,
                               [[0, 0, -13.5, -13.5, 13.5, 0],
                                [-13.5, 0, 0, 0, 13.5, 13.5]],
                               rtol=0, atol=1e-14)
    np.testing.assert_allclose(wrench.linear_wave,
                               [[0, 0, -22.5, -22.5, 22.5, 0]] * 2,
                               rtol=0, atol=1e-14)

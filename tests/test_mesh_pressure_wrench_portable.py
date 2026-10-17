"""Analytic force and moment checks for moving nonlinear-pressure facets."""

from pathlib import Path

import numpy as np

from wecsim.hydrodynamics.nonlinearHydro import (
    MeshPressures, mesh_pressure_wrenches, regular_wave_mesh_pressures,
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


def test_deep_water_facet_pressures_and_fk_correction_against_matlab():
    """Crest, trough, and heave must preserve dry-facet clipping."""
    source_file = (Path(__file__).parent / "fixtures/matlab"
                   / "deep_water_mesh_pressure/source.csv")
    source = np.loadtxt(source_file, delimiter=",", ndmin=2)
    assert source.shape == (3, 13)
    np.testing.assert_array_equal(source[:, :2], [[0, 0], [4, 0], [0, 1]])

    triangle = ((-1/3, -1/3), (2/3, -1/3), (-1/3, 2/3))
    vertices = np.array([[x, y, z] for z in (-2, -.25)
                         for x, y in triangle])
    faces = np.array([[0, 1, 2], [3, 4, 5]])
    poses = np.zeros((len(source), 6))
    poses[:, 2] = source[:, 1]
    pressure = regular_wave_mesh_pressures(
        vertices, faces, poses, source[:, 0],
        center_gravity=[0, 0, 0], rho=1000, gravity=9.81,
        water_depth=200, wave_number=.1, wave_height=1,
        wave_period=8, ramp_time=0, deep_water=True,
    )
    np.testing.assert_allclose(
        source[:, 2], .5 * np.cos(2 * np.pi * source[:, 0] / 8),
        rtol=0, atol=1e-15,
    )
    np.testing.assert_allclose(pressure.nonlinear_wave, source[:, 3:5],
                               rtol=0, atol=1e-8)
    np.testing.assert_allclose(pressure.linear_wave, source[:, 5:7],
                               rtol=0, atol=1e-8)

    # The shallow facet is wet at the crest, then dry at the trough and
    # after upward body motion. The mean-surface linear pressure remains.
    assert source[0, 4] > 0
    np.testing.assert_array_equal(source[1:, 4], [0, 0])
    np.testing.assert_allclose(pressure.hydrostatic,
                               [[19620, 2452.5], [19620, 0], [9810, 0]],
                               rtol=0, atol=1e-10)
    wrench = mesh_pressure_wrenches(vertices, faces, poses, pressure)
    np.testing.assert_allclose(
        wrench.nonlinear_wave - wrench.linear_wave, source[:, 7:13],
        rtol=0, atol=1e-8,
    )

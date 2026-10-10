"""Pair actual OSWEC flap pressures with the pinned MATLAB application."""

import os
from pathlib import Path

import numpy as np
import pytest
from scipy.io import loadmat
from scipy.optimize import brentq
import trimesh

from wecsim.bodyClass import BodyClass
from wecsim.nonlinearHydro import (
    mesh_pressure_wrenches,
    regular_wave_mesh_pressures,
)


def test_regular_mesh_pressure_uses_mean_and_moved_surfaces():
    vertices = np.array([[0, 0, -2], [1, 0, -2], [0, 1, -2]])
    faces = np.array([[0, 1, 2]])
    poses = np.array([[0, 0, 0, 0, 0, 0], [0, 0, 0.25, 0, 0, 0]])
    pressure = regular_wave_mesh_pressures(
        vertices, faces, poses, [0, 8], center_gravity=[0, 0, 0],
        rho=1000, gravity=9.81, water_depth=10, wave_number=0.1,
        wave_height=1, wave_period=8, ramp_time=0,
    )
    np.testing.assert_allclose(pressure.hydrostatic[:, 0], [19620, 17167.5])
    assert pressure.linear_wave[0, 0] == pressure.linear_wave[1, 0]
    assert pressure.nonlinear_wave[0, 0] != pressure.nonlinear_wave[1, 0]


@pytest.mark.skipif(
    not (os.getenv("WEC_SIM_MATLAB_OSWEC_PRESSURE_DIR")
         and os.getenv("WEC_SIM_OSWEC_FLAP_STL")),
    reason="pinned MATLAB OSWEC pressure output and flap STL not provided",
)
def test_published_oswec_nonlinear_visualization_pressures():
    reference = Path(os.environ["WEC_SIM_MATLAB_OSWEC_PRESSURE_DIR"])
    source = loadmat(reference / "source.mat")
    time = source["time"].ravel()
    body_time = source["bodyTime"].ravel()
    wave = source["wave"]
    poses = source["pose"]
    mesh = trimesh.load_mesh(os.environ["WEC_SIM_OSWEC_FLAP_STL"], process=False)
    assert len(mesh.faces) == 1042
    np.testing.assert_allclose(mesh.triangles_center, source["centers"],
                               rtol=0, atol=1e-9)
    np.testing.assert_allclose(time, body_time, rtol=0, atol=1e-12)
    np.testing.assert_allclose(time, wave[:, 0], rtol=0, atol=1e-12)
    assert len(time) == 1201 and time[-1] == 120

    depth = float(source["waterDepth"].item())
    source_k = float(source["waveNumber"].item())
    omega = 2 * np.pi / 8
    wave_number = brentq(
        lambda k: 9.81 * k * np.tanh(k * depth) - omega**2,
        1e-12, 1,
    )
    np.testing.assert_allclose(wave_number, source_k, rtol=0, atol=1e-10)
    assert int(source["deepWater"].item()) == 0
    assert depth == 10.9
    assert float(source["rho"].item()) == 1000
    assert float(source["gravity"].item()) == 9.81
    ramp = np.where(time >= 40, 1, (1 - np.cos(np.pi * time / 40)) / 2)
    elevation = 1.25 * ramp * np.cos(omega * time)
    np.testing.assert_allclose(elevation, wave[:, 1], rtol=0, atol=1e-10)

    pressure = regular_wave_mesh_pressures(
        mesh.vertices, mesh.faces, poses, time,
        center_gravity=source["cg"].ravel(), rho=1000, gravity=9.81,
        water_depth=depth, wave_number=wave_number,
        wave_height=2.5, wave_period=8, ramp_time=40,
    )
    for field, key in (
        (pressure.hydrostatic, "hydrostatic"),
        (pressure.nonlinear_wave, "nonlinearWave"),
        (pressure.linear_wave, "linearWave"),
    ):
        assert field.shape == (1201, 1042)
        np.testing.assert_allclose(field, source[key], rtol=0, atol=1e-8)



@pytest.mark.skipif(
    not (os.getenv("WEC_SIM_MATLAB_OSWEC_PRESSURE_DIR")
         and os.getenv("WEC_SIM_OSWEC_FLAP_STL")
         and os.getenv("WEC_SIM_OSWEC_H5")),
    reason="pinned MATLAB OSWEC forces, flap STL, and HDF5 not provided",
)
def test_published_oswec_nonlinear_visualization_forces():
    """Reconstruct logged forces from Python pressures on the source poses."""
    source = loadmat(Path(os.environ["WEC_SIM_MATLAB_OSWEC_PRESSURE_DIR"]) / "source.mat")
    mesh = trimesh.load_mesh(os.environ["WEC_SIM_OSWEC_FLAP_STL"], process=False)
    time = source["time"].ravel()
    poses = source["pose"]
    omega = 2 * np.pi / 8
    depth = float(source["waterDepth"].item())
    wave_number = brentq(
        lambda k: 9.81 * k * np.tanh(k * depth) - omega**2, 1e-12, 1,
    )
    pressure = regular_wave_mesh_pressures(
        mesh.vertices, mesh.faces, poses, time,
        center_gravity=source["cg"].ravel(), rho=1000, gravity=9.81,
        water_depth=depth, wave_number=wave_number,
        wave_height=2.5, wave_period=8, ramp_time=40,
    )
    ramp = np.where(time >= 40, 1, (1 - np.cos(np.pi * time / 40)) / 2)
    assert str(source["solver"].item()) == "ode4"
    np.testing.assert_allclose(mesh.face_normals, source["normals"],
                               rtol=0, atol=1e-12)
    np.testing.assert_allclose(mesh.area_faces, source["areas"].ravel(),
                               rtol=0, atol=1e-12)
    wrench = mesh_pressure_wrenches(mesh.vertices, mesh.faces, poses, pressure)
    hydrostatic = wrench.hydrostatic.copy()
    hydrostatic[:, 2] -= 127000 * 9.81
    np.testing.assert_allclose(source["forceRestoring"], -hydrostatic,
                               rtol=0, atol=2e-8)

    body = BodyClass(os.environ["WEC_SIM_OSWEC_H5"])
    body.bodyNumber = body.bodyTotal = 1
    body.readH5file()
    body.mass = 127000
    body.hydroStiffness = np.zeros((6, 6))
    body.viscDrag = {
        "Drag": np.zeros((6, 6)), "cd": np.zeros(6),
        "characteristicArea": np.zeros(6),
    }
    body.linearDamping = np.zeros((6, 6))
    body.hydroForcePre(
        omega, [0], 1, np.array([0.0]), [], 0.1, 1000, 9.81,
        "regular", np.vstack((time, np.zeros_like(time))),
        1, 1, 0, 2, 0,
    )
    re = np.asarray(body.hydroForce["fExt"]["re"])
    im = np.asarray(body.hydroForce["fExt"]["im"])
    linear_excitation = 1.25 * ramp[:, None] * (
        np.cos(omega * time)[:, None] * re
        - np.sin(omega * time)[:, None] * im
    )
    excitation = (linear_excitation
                  + ramp[:, None] * (wrench.nonlinear_wave - wrench.linear_wave))
    np.testing.assert_allclose(source["forceExcitation"], excitation,
                               rtol=0, atol=1e-6)
    radiation = source["velocity"] @ np.asarray(body.hydroForce["fDamping"]).T
    np.testing.assert_allclose(source["forceRadiationDamping"], radiation,
                               rtol=0, atol=1e-8)

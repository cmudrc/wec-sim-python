"""Physical invariants for the public instantaneous mesh-heave layout."""

from pathlib import Path

import numpy as np
import trimesh

from wecsim import RegularWave, WEC


HYDRO = Path(__file__).resolve().parents[1] / "examples/data/oswec.h5"


def test_submerged_mesh_drag_opposes_motion_without_a_wave(tmp_path):
    # The 1 m³ box stays fully submerged at the HDF5 body's -3.9 m CG.
    # Its equilibrium mass is its displaced water mass, so there is no
    # restoring force or wave correction to obscure the drag sign.
    geometry = tmp_path / "box.stl"
    trimesh.creation.box(extents=(1, 1, 1)).export(geometry)
    wec = WEC("submerged mesh heave")
    body = wec.body(
        "flap", HYDRO, inertia=(0, 1.85e6, 0),
        geometry_file=geometry, nonlinear_hydro="instantaneous",
        drag_coefficient=1, drag_area=1,
    )
    wec.coordinate("heave", body.move("heave"))

    def run(speed):
        return wec.run(
            RegularWave(0, 8), dt=.05, end_time=.5, ramp_time=0,
            initial_speed={"heave": speed}, rho=1000,
        )

    upward, downward = run(.2), run(-.2)
    for response in (upward, downward):
        speed = response.coordinates["heave"].velocity
        force = dict(response.raw.extra_outputs)
        np.testing.assert_allclose(
            force["body1_buoyancy_minus_weight"], 0, rtol=0, atol=1e-8,
        )
        np.testing.assert_array_equal(force["body1_nonlinear_fk_correction"], 0)
        np.testing.assert_allclose(
            force["body1_quadratic_drag_force"],
            -.5 * 1000 * speed * np.abs(speed), rtol=0, atol=1e-10,
        )
        assert np.all(force["body1_quadratic_drag_force"] * speed < 0)
        assert np.all(np.diff(np.abs(speed)) < 0)
        np.testing.assert_allclose(
            response.bodies["flap"].position[:, 2],
            -3.9 + response.coordinates["heave"].position,
            rtol=0, atol=1e-14,
        )

    np.testing.assert_allclose(
        downward.coordinates["heave"].position,
        -upward.coordinates["heave"].position, rtol=0, atol=1e-10,
    )
    np.testing.assert_allclose(
        downward.coordinates["heave"].velocity,
        -upward.coordinates["heave"].velocity, rtol=0, atol=1e-10,
    )

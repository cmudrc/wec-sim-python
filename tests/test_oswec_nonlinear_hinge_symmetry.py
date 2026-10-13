"""Portable rigid-hinge checks for the nonlinear OSWEC pressure solver."""

from pathlib import Path

import numpy as np
import trimesh

from wecsim.dynamics.hingePitch import solve_hinged_pitch_nonlinear_regular


HYDRO = Path(__file__).resolve().parents[1] / "examples/data/oswec.h5"


def test_symmetric_submerged_flap_has_no_unforced_drift(tmp_path):
    """Opposite PTO preloads must mirror a symmetric, wave-free hinge."""
    # A closed 1 m³ box at the fixture body's CG is fully below the surface.
    # Its 1,000 kg mass balances buoyancy, leaving no static pitch moment.
    mesh = trimesh.creation.box(extents=(1, 1, 1))
    mesh.apply_translation((0, 0, -3.9))
    geometry = tmp_path / "symmetric_flap.stl"
    mesh.export(geometry)

    def run(preload):
        return solve_hinged_pitch_nonlinear_regular(
            HYDRO, geometry, wave_height=0, wave_period=8,
            hinge_z=-8.9, body_mass=1000, pitch_inertia=1.85e6,
            pto_damping=1e5, pto_stiffness=1e8,
            pto_equilibrium=preload, dt=.05, end_time=1, ramp_time=0,
        )

    negative, neutral, positive = (run(preload) for preload in (-.1, 0, .1))
    assert np.max(np.abs(neutral.angle)) < 1e-12
    assert np.max(np.abs(neutral.angular_velocity)) < 1e-12
    assert positive.angle[-1] > .05
    assert negative.angle[-1] < -.05

    np.testing.assert_allclose(negative.angle, -positive.angle,
                               rtol=0, atol=1e-7)
    np.testing.assert_allclose(negative.angular_velocity,
                               -positive.angular_velocity, rtol=0, atol=1e-7)
    np.testing.assert_allclose(negative.pto_torque,
                               -positive.pto_torque, rtol=0, atol=2)
    np.testing.assert_allclose(negative.center_position[:, 0],
                               -positive.center_position[:, 0], rtol=0,
                               atol=1e-6)
    np.testing.assert_allclose(negative.center_position[:, 2],
                               positive.center_position[:, 2], rtol=0,
                               atol=1e-6)

    for response in (negative, neutral, positive):
        # The center follows a circle of radius 5 m about the fixed hinge.
        radius = np.hypot(response.center_position[:, 0],
                          response.center_position[:, 2] + 8.9)
        np.testing.assert_allclose(radius, 5, rtol=0, atol=1e-10)
        np.testing.assert_array_equal(response.excitation_force, 0)

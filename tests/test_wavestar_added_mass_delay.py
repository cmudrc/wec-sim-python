"""Audit applied added-mass feedback in the published WaveStar source prefix.

The saved force is the force applied inside the model. The body's reported
``forceAddedMass`` is reconstructed after the run with its original matrix.
"""

from pathlib import Path

import h5py
import numpy as np

from wecsim.dynamics.wavestar import WaveStarLinkage


ROOT = Path(__file__).parent / "fixtures/matlab"
FIXTURE = ROOT / "wavestar_added_mass_delay/source.npz"
PUBLISHED = ROOT / "wavestar_fixed_input_step/source.npz"
HYDRO = ROOT / "wavestar_fixed_input_step/wavestar.h5"


def test_source_added_mass_is_shifted_to_rigid_mass_and_delayed():
    with np.load(FIXTURE, allow_pickle=False) as saved, \
            np.load(PUBLISHED, allow_pickle=False) as published, \
            h5py.File(HYDRO) as hydro:
        time = saved["time"]
        pitch = saved["pitch"]
        speed = saved["speed"]
        acceleration = saved["acceleration"]
        applied = saved["applied_added_mass"]
        modified = saved["modified_added_mass"]
        original = saved["original_added_mass"]
        assert time.shape == pitch.shape == speed.shape == (200,)
        assert acceleration.shape == applied.shape == (200, 6)
        assert modified.shape == original.shape == (6, 6)
        np.testing.assert_allclose(time, np.arange(200) * .05,
                                   rtol=0, atol=1e-12)
        body = published["50ms_WECCCOMP_Nonlinear_Model_Predictive_body1"]
        np.testing.assert_allclose(pitch, body[:, 5], rtol=0, atol=1e-12)
        np.testing.assert_allclose(speed, body[:, 11], rtol=0, atol=1e-12)

        pinned_mass = np.asarray(
            hydro["body1/hydro_coeffs/added_mass/inf_freq"]
        ) * 1000
        np.testing.assert_allclose(original, pinned_mass, rtol=0, atol=1e-8)
        shifted_mass = 2 * np.trace(original[:3, :3])
        np.testing.assert_allclose(saved["adjusted_mass"],
                                   saved["rigid_mass"] + shifted_mass,
                                   rtol=0, atol=1e-9)
        np.testing.assert_allclose(modified[:3, :3],
                                   original[:3, :3] - shifted_mass * np.eye(3),
                                   rtol=0, atol=1e-8)
        np.testing.assert_allclose(np.diag(modified)[3:], 0,
                                   rtol=0, atol=1e-12)

        # The source's 100 ns transport delay prevents instantaneous
        # acceleration feedback. At saved 50 ms ticks, its applied force
        # follows a one-step history and then a linear extrapolation.
        np.testing.assert_allclose(applied[0], 0, rtol=0, atol=1e-12)
        np.testing.assert_allclose(applied[1], modified @ acceleration[0],
                                   rtol=0, atol=1e-10)
        extrapolated = (2 * acceleration[1:-1] - acceleration[:-2]) @ modified.T
        assert np.max(np.abs(applied[2:] - extrapolated)) < 1e-5


def test_mass_adjustment_changes_published_step_motion():
    with np.load(FIXTURE, allow_pickle=False) as saved:
        assert saved["mass_factor"] == 2
        assert saved["factor_one_pitch"].shape == (200,)
        assert saved["factor_one_speed"].shape == (200,)
        assert .3383 < saved["acceleration"][0, 4] < .3385
        assert .5191 < saved["factor_one_acceleration"][0, 4] < .5193
        pitch_change = np.max(np.abs(
            saved["pitch"] - saved["factor_one_pitch"]
        ))
        assert .0025 < pitch_change < .003


def test_bc_rod_and_motor_follow_published_linkage_kinematics():
    """Reconstruct both light body origins from the arm and B-to-C joints."""
    with np.load(FIXTURE, allow_pickle=False) as saved:
        angle = saved["pitch"]
        speed = saved["speed"]
        pitch_acceleration = saved["acceleration"][:, 4]
        linkage = WaveStarLinkage()
        base = np.asarray(linkage.pto_base)
        tip0 = np.asarray(linkage.pto_arm)
        tip = linkage.point_position(tip0, angle)
        tip_tangent = linkage.point_velocity(tip0, angle, 1.)
        rod0 = tip0 - base
        rod = tip - base
        theta = np.arctan2(
            rod0[1] * rod[:, 0] - rod0[0] * rod[:, 1],
            rod @ rod0,
        )
        rod_length_squared = np.sum(rod**2, axis=1)
        numerator = (rod[:, 1] * tip_tangent[:, 0]
                     - rod[:, 0] * tip_tangent[:, 1])
        theta_prime = numerator / rod_length_squared
        tip_second = -(tip - linkage.pivot)
        numerator_prime = (rod[:, 1] * tip_second[:, 0]
                           - rod[:, 0] * tip_second[:, 1])
        theta_second = (
            numerator_prime * rod_length_squared
            - 2 * numerator * np.sum(rod * tip_tangent, axis=1)
        ) / rod_length_squared**2

        def rotate(point):
            return np.stack((
                np.cos(theta) * point[0] + np.sin(theta) * point[1],
                -np.sin(theta) * point[0] + np.cos(theta) * point[1],
            ), axis=-1)

        rotated_tip = rotate(tip0)
        rotated_base = rotate(base)
        rod_position = tip - rotated_tip
        motor_position = base - rotated_base
        tip_rotation_tangent = np.stack((
            rotated_tip[:, 1], -rotated_tip[:, 0],
        ), axis=-1)
        base_rotation_tangent = np.stack((
            rotated_base[:, 1], -rotated_base[:, 0],
        ), axis=-1)
        rod_jacobian = tip_tangent - theta_prime[:, None] * tip_rotation_tangent
        motor_jacobian = -theta_prime[:, None] * base_rotation_tangent
        rod_second = (tip_second - theta_second[:, None] * tip_rotation_tangent
                      + theta_prime[:, None]**2 * rotated_tip)
        motor_second = (-theta_second[:, None] * base_rotation_tangent
                        + theta_prime[:, None]**2 * rotated_base)
        for number, position, jacobian, second in (
            (4, rod_position, rod_jacobian, rod_second),
            (5, motor_position, motor_jacobian, motor_second),
        ):
            np.testing.assert_allclose(saved[f"aux{number}_position"][:, [0, 2]],
                                       position, rtol=0, atol=1e-8)
            np.testing.assert_allclose(saved[f"aux{number}_position"][:, 4],
                                       theta, rtol=0, atol=1e-8)
            np.testing.assert_allclose(saved[f"aux{number}_velocity"][:, [0, 2]],
                                       jacobian * speed[:, None],
                                       rtol=0, atol=3e-8)
            np.testing.assert_allclose(saved[f"aux{number}_velocity"][:, 4],
                                       theta_prime * speed, rtol=0, atol=3e-8)
            np.testing.assert_allclose(
                saved[f"aux{number}_acceleration"][:, [0, 2]],
                jacobian * pitch_acceleration[:, None]
                + second * speed[:, None]**2,
                rtol=0, atol=2e-7,
            )
            np.testing.assert_allclose(
                saved[f"aux{number}_acceleration"][:, 4],
                theta_prime * pitch_acceleration + theta_second * speed**2,
                rtol=0, atol=2e-7,
            )

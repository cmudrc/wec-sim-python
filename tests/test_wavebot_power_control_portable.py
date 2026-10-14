"""Full published WaveBot controller trajectory against pinned MATLAB output."""

from pathlib import Path

import numpy as np
from scipy.io import loadmat

from wecsim.control.wavebotControl import run_wavebot_power_control
from wecsim.environment.irregularWave import pm_equal_energy_components


FIXTURE = Path(__file__).parent / "fixtures/matlab/wavebot_control"


def test_published_100s_adaptive_power_controller():
    """Regenerate the sea and all 13 gain updates without saved source states."""
    with np.load(FIXTURE / "published_100s.npz") as source:
        components = pm_equal_energy_components(
            FIXTURE / "waveBotBuoy.h5",
            significant_height=.254,
            peak_period=3.5,
            directions=np.array([0.]),
            spreading=np.array([1.]),
            phase=source["components"][:, 3, None],
        )
        np.testing.assert_allclose(
            np.column_stack((components.omega, components.spectral_amplitude,
                             components.d_omega)),
            source["components"][:, :3], rtol=0, atol=1e-12,
        )
        impedance = loadmat(FIXTURE / "Z_included.mat")["Z"]
        frequency = loadmat(FIXTURE / "f_vec.mat")["f_vec"].ravel()
        response = run_wavebot_power_control(
            FIXTURE / "waveBotBuoy.h5", components, impedance, frequency,
            source_linear_damping=True,
        )

        np.testing.assert_allclose(response.time, source["body_time"],
                                   rtol=0, atol=1e-8)
        np.testing.assert_allclose(response.time, source["wave"][:, 0],
                                   rtol=0, atol=1e-8)
        np.testing.assert_allclose(response.wave_elevation,
                                   source["wave"][:, 1], rtol=0, atol=1e-10)
        np.testing.assert_allclose(response.body_position[:, [0, 2, 4]],
                                   source["body_position"], rtol=0, atol=5e-4)
        np.testing.assert_allclose(response.body_velocity[:, [0, 2, 4]],
                                   source["body_velocity"], rtol=0, atol=2e-3)

        # Simulink records actuation on a separate grid; pair at body times.
        actuation = source["actuation"]
        command = np.column_stack([
            np.interp(response.time, actuation[:, 0], actuation[:, channel])
            for channel in (1, 2, 3)
        ])
        np.testing.assert_allclose(response.actuator_command, command,
                                   rtol=0, atol=2.)
        np.testing.assert_allclose(response.gains[0], 0, rtol=0, atol=1e-12)
        np.testing.assert_allclose(response.gains[1],
                                   [-1000, 0, -1000, 500, -200, 0],
                                   rtol=0, atol=1e-12)
        np.testing.assert_allclose(response.gains, source["gains"],
                                   rtol=0, atol=40.)
        np.testing.assert_allclose(response.excitation_estimates,
                                   source["excitation_estimates"],
                                   rtol=0, atol=.07)

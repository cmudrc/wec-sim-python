"""Independent WaveStar fault motion against the pinned MATLAB application."""

import os
from pathlib import Path

import numpy as np

from wecsim.environment.irregularWave import jonswap_equal_energy_components
from wecsim.control.wavestarFault import run_wavestar_fault_published


REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
LIVE_HYDRO = os.environ.get("WEC_SIM_WAVESTAR_H5")
LIVE = bool(REFERENCE and LIVE_HYDRO)
ROOT = Path(__file__).resolve().parent
PORTABLE = ROOT / "fixtures/matlab/wavestar_fault.npz"
PORTABLE_HYDRO = ROOT / "fixtures/matlab/wavestar_fault.h5"


def _source():
    if not LIVE:
        with np.load(PORTABLE, allow_pickle=False) as archive:
            return {name: archive[name] for name in archive.files}

    def csv(name):
        values = np.loadtxt(
            Path(REFERENCE) / f"WECCCOMP_FAULT_SOURCE_{name}.csv",
            delimiter=",", ndmin=2,
        )
        assert np.isfinite(values).all(), name
        return values

    float_body = csv("WECCCOMP_Fault_Implementation_body1")
    arm = csv("WECCCOMP_Fault_Implementation_body2")
    pto = csv("WECCCOMP_Fault_Implementation_pto1")
    true = csv("fault_true_position")
    measured = csv("fault_measured_position")
    assert true.shape == measured.shape == (141_201, 2)
    np.testing.assert_allclose(true[:, 0], np.arange(141_201) * .001,
                               rtol=0, atol=1e-8)
    np.testing.assert_allclose(measured[:, 0], true[:, 0],
                               rtol=0, atol=1e-12)
    return {
        "components": csv("components"),
        "sensor_noise": measured[:, 1] - true[:, 1],
        "sensor_dropout": measured[:, 1] == 0,
        "time": float_body[:, 0],
        "float_position": float_body[:, [1, 3, 5]],
        "float_velocity": float_body[:, [7, 9, 11]],
        "float_excitation": float_body[:, 19:25],
        "arm_position": arm[:, [1, 3, 5]],
        "arm_velocity": arm[:, [7, 9, 11]],
        "pto_stroke": pto[:, 3],
        "pto_speed": pto[:, 9],
        "pto_force": pto[:, 33],
        "radiation_force": csv("body1_forces")[:, 1:7],
        "wave_elevation": csv("wave")[:, 1],
    }


def test_independent_published_fault_trajectory():
    source = _source()
    components = source["components"]
    assert components.shape == (500, 4)
    assert source["sensor_noise"].shape == (141_201,)
    assert source["sensor_dropout"].shape == (141_201,)
    hydro = LIVE_HYDRO if LIVE else PORTABLE_HYDRO
    sea = jonswap_equal_energy_components(
        hydro, significant_height=.0625, peak_period=1.412,
        directions=np.array([0.]), spreading=np.array([1.]),
        phase=components[:, 3, None],
    )
    response = run_wavestar_fault_published(
        hydro, sea, source["sensor_noise"], source["sensor_dropout"],
    )
    assert response.time.shape == (14_121,)
    np.testing.assert_allclose(response.time, source["time"],
                               rtol=0, atol=1e-8)
    np.testing.assert_allclose(response.wave_elevation,
                               source["wave_elevation"],
                               rtol=0, atol=1e-10)
    np.testing.assert_allclose(response.excitation_force,
                               source["float_excitation"], rtol=0, atol=1e-7)
    np.testing.assert_allclose(response.angle,
                               source["float_position"][:, 2],
                               rtol=0, atol=2e-4)
    np.testing.assert_allclose(response.angular_speed,
                               source["float_velocity"][:, 2],
                               rtol=0, atol=.015)
    for position, velocity, expected_position, expected_velocity in (
        (response.float_position, response.float_velocity,
         source["float_position"], source["float_velocity"]),
        (response.arm_position, response.arm_velocity,
         source["arm_position"], source["arm_velocity"]),
    ):
        np.testing.assert_allclose(position[:, [0, 2]], expected_position[:, :2],
                                   rtol=0, atol=1e-4)
        np.testing.assert_allclose(velocity[:, [0, 2]], expected_velocity[:, :2],
                                   rtol=0, atol=.007)
        np.testing.assert_allclose(position[:, 4], expected_position[:, 2],
                                   rtol=0, atol=2e-4)
        np.testing.assert_allclose(velocity[:, 4], expected_velocity[:, 2],
                                   rtol=0, atol=.015)
    np.testing.assert_allclose(response.pto_stroke, source["pto_stroke"],
                               rtol=0, atol=4e-5)
    np.testing.assert_allclose(response.pto_speed, source["pto_speed"],
                               rtol=0, atol=.003)
    np.testing.assert_allclose(response.pto_force, source["pto_force"],
                               rtol=0, atol=.3)
    np.testing.assert_allclose(response.radiation_force,
                               source["radiation_force"], rtol=0, atol=.01)

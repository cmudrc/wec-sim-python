"""Reject damaged fitted radiation before advancing the WaveStar plant."""

from pathlib import Path
import shutil

import h5py
import numpy as np
import pytest

from wecsim.dynamics.wavestar import run_wavestar_published
from wecsim.environment.irregularWave import jonswap_equal_energy_components


ROOT = Path(__file__).parent / "fixtures/matlab"
HYDRO = ROOT / "wavestar_fixed_input_step/wavestar.h5"
SOURCE = ROOT / "wavestar_fixed_input_step/source.npz"


@pytest.mark.parametrize("damage", ["nonfinite_matrix", "fractional_order",
                                     "oversized_order"])
def test_wavestar_rejects_damaged_radiation_fit(tmp_path, damage):
    with np.load(SOURCE, allow_pickle=False) as source:
        sea = jonswap_equal_energy_components(
            HYDRO, significant_height=.1042, peak_period=1.836,
            directions=np.array([0.]), spreading=np.array([1.]),
            gamma=3.3, phase=source["50ms_components"][:, 3, None],
        )
    damaged = tmp_path / "wavestar.h5"
    shutil.copyfile(HYDRO, damaged)
    with h5py.File(damaged, "r+") as h5:
        fit = h5["body1/hydro_coeffs/radiation_damping/state_space"]
        if damage == "nonfinite_matrix":
            fit["A/all"][4, 4, 0, 0] = np.nan
        elif damage == "fractional_order":
            fit["it"][4, 4] = 2.5
        else:
            fit["it"][4, 4] = fit["A/all"].shape[2] + 1

    with pytest.raises(ValueError, match="finite six-DOF radiation fit"):
        run_wavestar_published(
            damaged, sea, dt=.05, end_time=.05, ramp_time=25,
            g=9.80665, solver="wecsim_ode8",
        )

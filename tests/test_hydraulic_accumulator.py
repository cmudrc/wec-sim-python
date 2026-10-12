"""Gas-volume pressure behavior for a hydraulic accumulator."""

import numpy as np
import pytest

from wecsim.components.hydraulic import GasChargedAccumulator


def test_inflow_compresses_gas_and_outflow_expands_it():
    accumulator = GasChargedAccumulator(8.5, 19_199_810.325)
    volumes = np.array([-.5, 0, .5])
    pressure = accumulator.pressure(volumes)
    assert pressure[0] < pressure[1] < pressure[2]
    np.testing.assert_allclose(pressure[1], accumulator.precharge_pressure)


def test_accumulator_rejects_exhausted_gas_volume():
    accumulator = GasChargedAccumulator(8.5, 19_199_810.325)
    with pytest.raises(ValueError, match="gas volume"):
        accumulator.pressure(8.5)
    with pytest.raises(ValueError, match="positive"):
        GasChargedAccumulator(0, 19_199_810.325)

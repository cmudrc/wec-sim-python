"""Validate switch clocks at the direct variable-heave boundary."""

import numpy as np
import pytest

from wecsim.hydrodynamics.variableHydro import (
    HeaveHydroState,
    integrate_variable_heave,
)


@pytest.mark.parametrize("switch", [True, np.bool_(True), "0.1"])
def test_direct_heave_rejects_switch_values_coerced_to_numbers(switch):
    state = HeaveHydroState(
        mass=1000, displaced_volume=1, equilibrium_z=0,
        added_mass=100, radiation_damping=10, hydrostatic_stiffness=100,
        excitation_real=0, excitation_imaginary=0,
    )
    with pytest.raises(ValueError, match="switch times must be real numbers"):
        integrate_variable_heave(
            (state, state), (switch,), dt=.1, end_time=.2,
            height=0, period=8, ramp_time=0, rho=1000, g=9.81,
            pto_stiffness=0, pto_damping=0,
        )

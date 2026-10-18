"""Validate switch clocks at the direct variable-heave boundary."""

from pathlib import Path

import numpy as np
import pytest

from wecsim import HydroState, RegularWave, WEC, WorldPoint
from wecsim.dynamics.caseDynamics import run_case
from wecsim.hydrodynamics.variableHydro import (
    HeaveHydroState,
    integrate_variable_heave,
)


SPHERE = Path(__file__).resolve().parents[1] / "examples/data/sphere.h5"


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


def test_raw_case_rejects_boolean_switch_time():
    wec = WEC("raw variable hydro")
    sphere = wec.variable_body(
        "sphere", (HydroState(SPHERE), HydroState(SPHERE)),
        switch_times=(.1,),
    )
    wec.coordinate("heave", sphere.move("heave"))
    wec.pto("damper", WorldPoint(0, 0, 2), sphere.at(0, 0, 0),
            axis=(0, 0, 1), damping=200_000)
    case = wec.to_case(RegularWave(0, 8), dt=.1, end_time=.2, ramp_time=0)
    case["bodies"][0]["variable_hydro"]["switch_times"] = [True]
    with pytest.raises(ValueError, match="switch times must be real numbers"):
        run_case(case)

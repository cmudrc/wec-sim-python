"""Pair a fixed monopile with current against pinned MATLAB WEC-Sim."""

from hashlib import sha256
import os
from pathlib import Path

import numpy as np
import pytest

from wecsim import Current, PMWave, WEC
from wecsim.environment.irregularWave import pm_equal_energy_components
from wecsim.hydrodynamics.morison import finite_depth_wavenumber


REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
PORTABLE = Path(__file__).parent / "fixtures/matlab/fixed_morison_current/source.npz"
PORTABLE_SHA256 = "203e0e7e702655b18c98b47bbc32ea2cab4f447c4859a0f95a34b4f75d972c07"


def _source_arrays(profile):
    if REFERENCE:
        source = Path(REFERENCE)
        prefix = f"MORISON_FIXED_CURRENT_{profile}"
        load = lambda suffix: np.loadtxt(source / f"{prefix}_{suffix}.csv",
                                        delimiter=",")
        return (load("components"), load("wave"), load("morison_force"),
                load("body1"), load("body2"))
    assert sha256(PORTABLE.read_bytes()).hexdigest() == PORTABLE_SHA256
    with np.load(PORTABLE, allow_pickle=False) as saved:
        return (saved["components"].copy(), saved["wave"].copy(),
                saved[f"{profile}_force"].copy(),
                saved[f"{profile}_body1"].copy(), saved["body2"].copy())


def test_public_fixed_morison_current_reaches_force_solver():
    wec = WEC("fixed current smoke test")
    body = wec.fixed_body(
        "pile", center_gravity=(0, 0, -5), mass="equilibrium",
        inertia=(1, 1, 1), volume=1,
    )
    wec.morison_element(
        body, point=body.at(0, 0, 0), drag_coefficient=(1, 1, 1),
        added_mass_coefficient=(0, 0, 0), area=(2, 2, 2), volume=1,
    )
    settings = dict(seed=1, directions=(30,), spreading=(1,),
                    frequency_count=4, frequency_range=(.1, 2),
                    water_depth=30)
    still = wec.run(PMWave(2, 5, **settings), dt=.1, end_time=.2,
                    ramp_time=.1)
    with_current = wec.run(
        PMWave(2, 5, **settings, current=Current(.8, 45, "power", 30)),
        dt=.1, end_time=.2, ramp_time=.1,
    )
    np.testing.assert_array_equal(with_current.wave_elevation,
                                  still.wave_elevation)
    assert np.max(np.abs(with_current.body_forces["pile"]
                             - still.body_forces["pile"])) > 1


def _check(actual, expected, tolerance, label):
    actual, expected = np.asarray(actual), np.asarray(expected)
    assert actual.shape == expected.shape, label
    assert np.isfinite(actual).all() and np.isfinite(expected).all(), label
    error = np.max(np.abs(actual - expected))
    assert error < tolerance, f"{label}: {error:.6g} exceeds {tolerance}"


@pytest.mark.parametrize("profile", ["uniform", "power", "linear"])
def test_fixed_monopile_current_profiles_against_matlab(tmp_path, profile):
    components, wave, force, monopile, tower = _source_arrays(profile)
    assert components.shape == (500, 5)
    assert wave.shape == (2001, 2)
    assert force.shape == (2001, 7)
    assert monopile.shape == tower.shape == (2001, 25)
    _check(monopile[:, 13:19], -force[:, 1:7], 1e-8,
           "MATLAB total and Morison force sign")

    phases = tmp_path / f"{profile}_phases.csv"
    np.savetxt(phases, components[:, 4:5], delimiter=",")
    sea = pm_equal_energy_components(
        None, significant_height=2, peak_period=5,
        directions=(30,), spreading=(1,), frequency_range=(.001, 10),
        phase=components[:, 4:5],
    )
    _check(sea.omega, components[:, 0], 1e-12, "frequency")
    _check(sea.spectral_amplitude, components[:, 1], 1e-12, "PM spectrum")
    _check(sea.d_omega, components[:, 2], 1e-12, "frequency width")
    _check(finite_depth_wavenumber(sea.omega, water_depth=30),
           components[:, 3], 1e-12, "wavenumber")

    wec = WEC("fixed monopile with current")
    body = wec.fixed_body(
        "monopile", center_gravity=(0, 0, -15),
        mass="equilibrium", inertia=(1.25e9, 1.25e9, .15e9),
        volume=np.pi * 10**2 * 30,
    )
    wec.fixed_body(
        "tower", center_gravity=(0, 0, 25), mass=1031930,
        inertia=(9.66e8, 9.66e8, .132e8),
    )
    wec.morison_element(
        body, point=body.at(0, 0, 10),
        drag_coefficient=(1, 1, 1),
        added_mass_coefficient=(1, 1, 1),
        area=(300, 300, np.pi * 10**2 / 4),
        volume=np.pi * 10**2 * 30,
    )
    result = wec.run(
        PMWave(2, 5, phase_file=phases, directions=(30,), spreading=(1,),
               frequency_range=(.001, 10), water_depth=30,
               current=Current(.8, 45, profile, 30)),
        dt=.01, end_time=20, ramp_time=10, rho=1025,
    )
    _check(result.time, wave[:, 0], 1e-12, "time")
    _check(result.wave_elevation, wave[:, 1], 1e-10, "wave elevation")
    _check(result.body_forces["monopile"], -force[:, 1:7], 1e-3,
           f"six-component {profile} current Morison force")
    _check(result.body_forces["tower"], np.zeros((2001, 6)), 1e-12,
           "tower force")
    _check(result.bodies["monopile"].position, monopile[:, 1:7], 1e-12,
           "fixed monopile position")
    _check(result.bodies["tower"].position, tower[:, 1:7], 1e-12,
           "fixed tower position")
    _check(result.bodies["monopile"].velocity, monopile[:, 7:13], 1e-12,
           "fixed monopile velocity")
    _check(result.bodies["tower"].velocity, tower[:, 7:13], 1e-12,
           "fixed tower velocity")

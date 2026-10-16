"""Reject case inputs that would silently change RM3 physics or geometry."""

from pathlib import Path

import numpy as np
import pytest

from wecsim import NoWave, PMWave, RegularWave, WEC
from wecsim.dynamics.caseDynamics import run_case


ROOT = Path(__file__).resolve().parents[1]
RM3 = ROOT / "examples/data/rm3.h5"
SPHERE = ROOT / "examples/data/sphere.h5"
OSWEC = ROOT / "examples/data/oswec.h5"


def _case(wave_mode):
    wec = WEC("RM3 physical boundaries")
    float_body = wec.body("float", RM3, inertia=(0, 21_306_090.66, 0))
    spar = wec.body("spar", RM3, inertia=(0, 94_407_091.24, 0))
    wec.floating_joint(float_body, spar, damping=1_200_000)
    if wave_mode == "regular":
        wave = RegularWave(2.5, 8)
    elif wave_mode == "pm":
        wave = PMWave(2.5, 8, seed=1, frequency_count=8)
    else:
        wave = NoWave()
    return wec.to_case(wave, dt=.1, end_time=.2,
                       ramp_time=None if wave_mode == "none" else 0)


def _replace(case, path, value):
    target = case
    for key in path[:-1]:
        target = target[key]
    target[path[-1]] = value


@pytest.mark.parametrize("wave_mode, path, value, message", (
    ("regular", ("simulation", "dt"), True, "simulation.dt must be numeric"),
    ("regular", ("simulation", "rho"), -1, "simulation.rho is outside"),
    ("regular", ("simulation", "end_time"), "bad", "simulation.end_time must be numeric"),
    ("regular", ("body_to_body",), "yes", "body_to_body must be a boolean"),
    ("regular", ("bodies",), [], "bodies must be a nonempty list"),
    ("regular", ("bodies", 1, "hydro_file"), str(SPHERE), "one shared HDF5"),
    ("regular", ("bodies", 1, "inertia"), [0, 1, 0], "mass and pitch_inertia"),
    ("regular", ("bodies", 1, "mass"), 1000, "require equilibrium mass"),
    ("regular", ("constraint", "coordinates"), [], "unsupported floating-joint"),
    ("regular", ("pto", "kind"), "heave", "relative_heave PTO"),
    ("regular", ("pto", "damping"), -1, "pto.damping is outside"),
    ("regular", ("pto",), {"kind": "relative_heave", "damping": 1e6,
                             "stiffness": 0, "pretension": 100},
     "pretension needs positive stiffness"),
    ("regular", ("mooring",), {"kind": "joint_surge_spring",
                                 "stiffness": 100, "session": object()},
     "joint surge spring needs only stiffness"),
    ("regular", ("mooring",), {"kind": "moor_dyn", "session": object()},
     "MoorDyn needs a native session"),
    ("regular", ("mooring",), {"kind": "unsupported"},
     "unsupported floating-joint mooring"),
    ("regular", ("wave", "direction"), 30, "support 0-degree waves"),
    ("regular", ("wave", "seed"), 1, "regular waves use height"),
    ("regular", ("simulation", "radiation_memory"), 1,
     "regular-wave floating-joint dynamics use constant radiation"),
    ("regular", ("simulation", "radiation_method"), "convolution",
     "regular waves need constant radiation"),
    ("pm", ("wave", "directions"), [30.], "one 0-degree wave heading"),
    ("pm", ("simulation", "radiation_method"), "constant",
     "radiation memory needs convolution"),
    ("none", ("simulation", "ramp_time"), 1,
     "no-wave floating joint has no wave or ramp settings"),
    ("none", ("wave", "height"), 1,
     "no-wave floating joint has no wave or ramp settings"),
))
def test_rm3_rejects_inconsistent_case_inputs(wave_mode, path, value, message):
    case = _case(wave_mode)
    _replace(case, path, value)
    with pytest.raises(ValueError, match=message):
        run_case(case, base_dir=ROOT)


def test_unmodified_rm3_case_still_runs_and_dissipates_pto_work():
    response = run_case(_case("regular"), base_dir=ROOT)
    assert response.time.shape == (3,)
    assert np.isfinite(response.body_position).all()
    np.testing.assert_allclose(
        response.pto_absorbed_power,
        1_200_000 * response.pto_velocity**2,
        rtol=0, atol=1e-8,
    )


def _hinge_case():
    wec = WEC("OSWEC physical boundaries")
    flap = wec.body("flap", OSWEC, mass=127_000,
                    inertia=(0, 1.85e6, 0))
    base = wec.fixed_body("base", center_gravity=(0, 0, -10.9), mass=999)
    wec.fixed_hinge(flap, base, damping=12_000, stiffness=500)
    return wec.to_case(RegularWave(.1, 8), dt=.1, end_time=.2, ramp_time=0)


@pytest.mark.parametrize("path, value, message", (
    (("wave", "type"), "none", "no-wave fixed hinge accepts no wave settings"),
    (("bodies", 0, "inertia"), [0, 1.85e6, 0],
     "mass, pitch inertia, and optional mesh hydro"),
    (("bodies", 0, "nonlinear_hydro"), "instantaneous",
     "nonlinear fixed hinge needs instantaneous hydro and a mesh"),
    (("constraint", "initial_coordinate"), {"surge": 1},
     "unknown coordinates"),
    (("pto",), {"kind": "pitch", "damping": 12_000},
     "fixed nonhydrodynamic base needs pto.location"),
    (("simulation", "added_mass_scheme"), "simulink_delay",
     "regular fixed-hinge dynamics need implicit added mass"),
    (("simulation", "solver"), "wecsim_ode4",
     "regular source ode4 needs nonlinear mesh pressures"),
    (("wave", "seed"), 1, "regular waves use height"),
    (("wave", "direction"), 30,
     "fixed-hinge regular dynamics currently support 0-degree waves"),
    (("simulation", "radiation_memory"), 1,
     "regular fixed-hinge dynamics use constant-frequency radiation"),
))
def test_oswec_rejects_inconsistent_hinge_inputs(path, value, message):
    case = _hinge_case()
    _replace(case, path, value)
    with pytest.raises(ValueError, match=message):
        run_case(case, base_dir=ROOT)


def test_unmodified_oswec_case_still_runs_with_its_torsional_pto():
    response = run_case(_hinge_case(), base_dir=ROOT)
    assert response.time.shape == (3,)
    np.testing.assert_allclose(
        response.pto_force,
        -12_000 * response.body_velocity[:, 0, 4]
        - 500 * response.body_position[:, 0, 4],
        rtol=0, atol=1e-9,
    )

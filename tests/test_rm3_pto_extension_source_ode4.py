"""Gate both published RM3 free decays through public source-style ode4."""

from pathlib import Path

import numpy as np
import pytest

from wecsim import NoWave, RegularWave, WEC


ROOT = Path(__file__).resolve().parents[1]
HYDRO = ROOT / "examples/data/rm3_end_stops.h5"
MATLAB = ROOT / "tests/fixtures/matlab/rm3_pto_extension_step.npz"


def _wec(*, added_mass_scheme="simulink_delay", damping=0):
    wec = WEC("RM3 PTO Extension")
    float_body = wec.body("float", HYDRO,
                          inertia=(0, 21_306_090.66, 0))
    spar = wec.body("spar", HYDRO, inertia=(0, 94_407_091.24, 0))
    wec.floating_joint(float_body, spar, damping=damping,
                       added_mass_scheme=added_mass_scheme)
    return wec


@pytest.mark.parametrize("name,active,column", (
    ("float", "float", 1), ("spar", "spar", 3),
))
@pytest.mark.parametrize("dt,suffix", (
    (.1, "0p1"), (.05, "0p05"), (.025, "0p025"),
    (.0125, "0p0125"),
))
def test_full_duration_published_source_ode4_motion(name, active, column,
                                                    dt, suffix):
    initial = {"float_heave": 5} if name == "float" else {"spar_heave": -5}
    result = _wec().run(
        NoWave(), dt=dt, end_time=30, radiation_memory=60,
        initial_coordinate=initial, solver="wecsim_ode4",
    )
    with np.load(MATLAB, allow_pickle=False) as archived:
        source = archived[f"{name}_{suffix}"].copy()
    assert source.shape == (round(30 / dt) + 1, 9)
    np.testing.assert_allclose(result.time, source[:, 0], rtol=0, atol=1e-10)
    body = result.bodies[active]
    # Each source trajectory is advanced from HDF5 and initial extension;
    # MATLAB positions, accelerations, and forces do not enter this solve.
    assert np.max(np.abs(body.position[:, 2] - source[:, column])) < 1e-8
    assert np.max(np.abs(body.velocity[:, 2] - source[:, column + 1])) < 1e-8
    inactive = result.bodies["spar" if name == "float" else "float"]
    inactive_column = 3 if name == "float" else 1
    assert np.max(np.abs(inactive.position[:, 2]
                         - source[:, inactive_column])) < 1e-8
    assert np.max(np.abs(inactive.velocity[:, 2]
                         - source[:, inactive_column + 1])) < 1e-8
    pto = result.ptos["relative_heave"]
    assert np.max(np.abs(pto.stroke - source[:, 5])) < 1e-8
    assert np.max(np.abs(pto.velocity - source[:, 6])) < 1e-8
    assert np.max(np.abs(pto.force - source[:, 7])) < 1e-7
    assert np.max(np.abs(pto.absorbed_power - source[:, 8])) < 1e-7
    assert result.wave_elevation is None


def test_source_ode4_is_explicit_and_limited_to_the_published_layout():
    initial = {"float_heave": 5}
    physical = _wec(added_mass_scheme="implicit").run(
        NoWave(), dt=.1, end_time=30, radiation_memory=60,
        initial_coordinate=initial,
    )
    with np.load(MATLAB, allow_pickle=False) as archived:
        published = archived["float_0p1"][:, 1].copy()
    assert np.max(np.abs(physical.bodies["float"].position[:, 2]
                         - published)) > .08
    for wec, wave in (
        (_wec(added_mass_scheme="implicit"), NoWave()),
        (_wec(damping=1), NoWave()),
        (_wec(), RegularWave(2.5, 8)),
    ):
        with pytest.raises(ValueError, match="source CIC ode4"):
            wec.run(wave, dt=.1, end_time=1, radiation_memory=60,
                    initial_coordinate=initial, solver="wecsim_ode4")

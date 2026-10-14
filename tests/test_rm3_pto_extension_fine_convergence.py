"""Bound the two published RM3 PTO-extension decays with finer source steps."""

from pathlib import Path

import numpy as np
import pytest

from wecsim import NoWave, WEC


ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests/fixtures/matlab/rm3_pto_extension_00078125/source.npz"
HYDRO = ROOT / "examples/data/rm3_end_stops.h5"


@pytest.mark.parametrize("name,column,source_limit,python_limit", (
    ("float", 1, (.001, .003), (.0002, .001)),
    ("spar", 3, (.0003, .0002), (.00002, .0002)),
))
def test_refined_source_and_independent_physical_motion(
    name, column, source_limit, python_limit,
):
    """Use the same full-duration screens declared before the local run."""
    with (np.load(ROOT / "tests/fixtures/matlab/rm3_pto_extension_00625.npz") as old,
          np.load(ROOT / "tests/fixtures/matlab/rm3_pto_extension_003125.npz") as prior,
          np.load(FIXTURE) as latest):
        c = old[f"{name}_0p00625"]
        d = prior[f"{name}_0p003125"]
        e = latest[f"{name}_0p0015625"]
        f = latest[f"{name}_0p00078125"]
        for dt, values in ((.00625, c), (.003125, d),
                           (.0015625, e), (.00078125, f)):
            assert values.shape == (round(30 / dt) + 1, 9)
            np.testing.assert_allclose(values[:, 0],
                                       np.arange(len(values)) * dt,
                                       rtol=0, atol=1e-9)
            assert np.isfinite(values).all()
            assert np.max(np.abs(values[:, 7:9])) < 1e-7

        channel = slice(column, column + 2)
        old_change = np.max(np.abs(c[:, channel] - d[::2, channel]), axis=0)
        mid_change = np.max(np.abs(d[:, channel] - e[::2, channel]), axis=0)
        final_change = np.max(np.abs(e[:, channel] - f[::2, channel]), axis=0)
        assert np.all(mid_change < old_change)
        assert np.all(final_change < mid_change)
        assert np.all(final_change < source_limit), (name, final_change)

        # The oldest source history is sampled at 0.00625 s. Compare all
        # first-order extrapolations on that common grid, without interpolation.
        old_extrap = 2 * d[::2, channel] - c[:, channel]
        mid_extrap = 2 * e[::4, channel] - d[::2, channel]
        final_extrap = 2 * f[::8, channel] - e[::4, channel]
        old_spread = np.max(np.abs(mid_extrap - old_extrap), axis=0)
        final_spread = np.max(np.abs(final_extrap - mid_extrap), axis=0)
        assert np.all(final_spread < old_spread), (name, final_spread)

        wec = WEC(f"RM3 PTO extension {name} convergence diagnostic")
        float_body = wec.body("float", HYDRO,
                              inertia=(0, 21_306_090.66, 0))
        spar = wec.body("spar", HYDRO,
                        inertia=(0, 94_407_091.24, 0))
        wec.floating_joint(float_body, spar, added_mass_scheme="implicit")
        initial = {"float_heave": 5} if name == "float" else {"spar_heave": -5}
        result = wec.run(NoWave(), dt=.00078125, end_time=30,
                         radiation_memory=60, initial_coordinate=initial)
        body = result.bodies[name]
        python = np.column_stack((body.position[::8, 2],
                                  body.velocity[::8, 2]))
        gap = np.max(np.abs(python - final_extrap), axis=0)
        assert np.all(gap < final_spread), (name, gap, final_spread)
        assert np.all(gap < python_limit), (name, gap, python_limit)

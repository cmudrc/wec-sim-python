"""Separate RM3 PTO-extension source step sensitivity from Python motion."""

from pathlib import Path

import numpy as np
import pytest

from wecsim import NoWave, WEC


ROOT = Path(__file__).resolve().parents[1]
# This pinned common-RM3 BEMIO file has the same active convolution,
# added-mass, restoring, volume, and CG data as the PTO-extension baseline.
HYDRO = ROOT / "examples/data/rm3_end_stops.h5"
MATLAB = ROOT / "tests/fixtures/matlab/rm3_pto_extension_step.npz"
MATLAB_FINER = ROOT / "tests/fixtures/matlab/rm3_pto_extension_00625.npz"
MATLAB_FINEST = ROOT / "tests/fixtures/matlab/rm3_pto_extension_003125.npz"


def _run(name, dt, scheme):
    wec = WEC(f"RM3 PTO extension {name}")
    float_body = wec.body("float", HYDRO,
                          inertia=(0, 21_306_090.66, 0))
    spar = wec.body("spar", HYDRO, inertia=(0, 94_407_091.24, 0))
    wec.floating_joint(float_body, spar, added_mass_scheme=scheme)
    initial = {"float_heave": 5} if name == "float" else {"spar_heave": -5}
    return wec.run(NoWave(), dt=dt, end_time=30, radiation_memory=60,
                   initial_coordinate=initial)


@pytest.mark.parametrize("name,column,position_limit,speed_limit", (
    ("float", 1, .012, .03),
    ("spar", 3, .003, .001),
))
def test_physical_pto_extension_step_refines_faster_than_matlab(
    name, column, position_limit, speed_limit,
):
    """Bound the physical path, then pair the opt-in source convention."""
    with np.load(MATLAB, allow_pickle=False) as source:
        samples = {suffix: source[f"{name}_{suffix}"]
                   for suffix in ("0p1", "0p05", "0p025", "0p0125")}
        for suffix, dt in (("0p1", .1), ("0p05", .05),
                           ("0p025", .025), ("0p0125", .0125)):
            assert samples[suffix].shape == (round(30 / dt) + 1, 9)
            np.testing.assert_allclose(samples[suffix][:, 0],
                                       np.arange(len(samples[suffix])) * dt,
                                       rtol=0, atol=1e-10)
        source_step = np.max(np.abs(
            samples["0p025"][:, column:column + 2]
            - samples["0p0125"][::2, column:column + 2],
        ), axis=0)
        prior_source_step = np.max(np.abs(
            samples["0p05"][:, column:column + 2]
            - samples["0p025"][::2, column:column + 2],
        ), axis=0)
        first_source_step = np.max(np.abs(
            samples["0p1"][:, column:column + 2]
            - samples["0p05"][::2, column:column + 2],
        ), axis=0)
        assert np.all(prior_source_step < first_source_step)
        assert np.all(source_step > .45 * prior_source_step)

        coarse = _run(name, .025, "implicit")
        fine = _run(name, .0125, "implicit")
        body = fine.bodies[name]
        coarse_body = coarse.bodies[name]
        physical_step = np.array([
            np.max(np.abs(coarse_body.position[:, 2]
                          - body.position[::2, 2])),
            np.max(np.abs(coarse_body.velocity[:, 2]
                          - body.velocity[::2, 2])),
        ])
        assert np.all(source_step > 5 * physical_step), (
            name, source_step, physical_step,
        )
        physical_gap = np.array([
            np.max(np.abs(body.position[:, 2]
                          - samples["0p0125"][:, column])),
            np.max(np.abs(body.velocity[:, 2]
                          - samples["0p0125"][:, column + 1])),
        ])
        print(f"{name}: physical self-change {physical_step}; "
              f"MATLAB self-change {source_step}; "
              f"physical/source gap {physical_gap}")
        assert physical_gap[0] < position_limit
        assert physical_gap[1] < speed_limit
        assert fine.wave_elevation is None
        assert np.max(np.abs(fine.ptos["relative_heave"].force)) < 1e-8

        source_scheme = _run(name, .0125, "simulink_delay")
        source_body = source_scheme.bodies[name]
        source_gap = np.array([
            np.max(np.abs(source_body.position[:, 2]
                          - samples["0p0125"][:, column])),
            np.max(np.abs(source_body.velocity[:, 2]
                          - samples["0p0125"][:, column + 1])),
        ])
        print(f"{name}: opt-in source-delay gap {source_gap}")
        limits = (.009, .017) if name == "float" else (.0015, .0004)
        assert np.all(source_gap < limits), (name, source_gap)


@pytest.mark.parametrize("name,column,position_limit,speed_limit", (
    ("float", 1, .012, .03),
    ("spar", 3, .003, .001),
))
def test_local_matlab_next_step_and_independent_physical_motion(
    name, column, position_limit, speed_limit,
):
    """Keep the new 30 s source refinement and physical pair in regular CI.

    The source comparison tests numerical stability, not the published 0.1 s
    trajectory. The float remains above the predeclared 1 mm / 3 mm/s screen.
    """
    with (np.load(MATLAB, allow_pickle=False) as archived,
          np.load(MATLAB_FINER, allow_pickle=False) as local):
        previous = archived[f"{name}_0p0125"]
        older = archived[f"{name}_0p025"]
        finer = local[f"{name}_0p00625"]
        assert previous.shape == (2401, 9)
        assert finer.shape == (4801, 9)
        np.testing.assert_allclose(finer[::2, 0], previous[:, 0],
                                   rtol=0, atol=1e-12)
        assert np.isfinite(finer).all()
        prior_change = np.max(np.abs(
            previous[::2, column:column + 2]
            - older[:, column:column + 2],
        ), axis=0)
        new_change = np.max(np.abs(
            finer[::2, column:column + 2]
            - previous[:, column:column + 2],
        ), axis=0)
        assert np.all(new_change < prior_change), (
            name, prior_change, new_change,
        )
        if name == "float":
            assert new_change[0] > .001 and new_change[1] > .003

        physical = _run(name, .00625, "implicit")
        body = physical.bodies[name]
        gap = np.max(np.abs(np.column_stack((
            body.position[:, 2], body.velocity[:, 2],
        )) - finer[:, column:column + 2]), axis=0)
        print(f"{name}: source step change {new_change}; "
              f"independent Python/source gap {gap}")
        assert gap[0] < position_limit
        assert gap[1] < speed_limit


@pytest.mark.parametrize("name,column", (("float", 1), ("spar", 3)))
def test_local_matlab_extrapolation_trend_and_physical_diagnostic(name, column):
    """Track the declared source convergence trend, not published-step parity."""
    with (np.load(MATLAB, allow_pickle=False) as archived,
          np.load(MATLAB_FINER, allow_pickle=False) as finer,
          np.load(MATLAB_FINEST, allow_pickle=False) as finest):
        a = archived[f"{name}_0p025"]
        b = archived[f"{name}_0p0125"]
        c = finer[f"{name}_0p00625"]
        d = finest[f"{name}_0p003125"]
        assert d.shape == (9601, 9)
        np.testing.assert_allclose(d[:, 0], np.arange(9601) * .003125,
                                   rtol=0, atol=1e-10)
        assert np.isfinite(d).all()
        assert np.max(np.abs(d[:, 7:9])) < 1e-7  # Published PTO is unforced.

        previous_step = np.max(np.abs(
            b[:, column:column + 2] - c[::2, column:column + 2]), axis=0)
        new_step = np.max(np.abs(
            c[:, column:column + 2] - d[::2, column:column + 2]), axis=0)
        assert np.all(new_step < previous_step), (name, previous_step, new_step)

        # Formal first-order extrapolations share the 0.025 s output grid.
        first = 2 * b[::2, column:column + 2] - a[:, column:column + 2]
        second = 2 * c[::4, column:column + 2] - b[::2, column:column + 2]
        third = 2 * d[::8, column:column + 2] - c[::4, column:column + 2]
        previous_inconsistency = np.max(np.abs(second - first), axis=0)
        new_inconsistency = np.max(np.abs(third - second), axis=0)
        assert np.all(new_inconsistency < previous_inconsistency), (
            name, previous_inconsistency, new_inconsistency)

        physical = _run(name, .003125, "implicit").bodies[name]
        python = np.column_stack((physical.position[::8, 2],
                                  physical.velocity[::8, 2]))
        physical_gap = np.max(np.abs(python - third), axis=0)
        # The independent run lies inside the remaining source-extrapolation
        # spread; this does not establish that the source has converged.
        assert np.all(physical_gap < new_inconsistency), (
            name, physical_gap, new_inconsistency)
        if name == "float":
            assert new_step[0] > .001 and new_step[1] > .003
        print(f"{name}: source step {new_step}; extrapolation spread "
              f"{new_inconsistency}; Python/extrapolation {physical_gap}")

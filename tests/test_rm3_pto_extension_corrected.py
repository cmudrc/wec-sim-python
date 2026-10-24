"""Gate the adopted RM3 PTO-extension decays against a corrected source."""

import os
from pathlib import Path

import h5py
import numpy as np
import pytest

from wecsim import NoWave, WEC


ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "tests/fixtures/matlab"
CORRECTED = FIXTURES / "rm3_pto_extension_corrected/source.npz"
BODY_AXES = FIXTURES / "rm3_pto_extension_corrected/body_axes.npz"
PRIOR = FIXTURES / "rm3_pto_extension_00078125/source.npz"
HYDRO = ROOT / "examples/data/rm3_end_stops.h5"
FRESH = os.environ.get("WEC_SIM_RM3_PTO_CORRECTED_DIR")
FRESH_HYDRO = os.environ.get("WEC_SIM_RM3_PTO_H5")
SELECTED = os.environ.get("WEC_SIM_RM3_PTO_CASE")

# The source_limit values were predeclared for MATLAB source self-change
# before the final two source runs. This corrected-baseline adoption also
# uses them as direct physical Python/source acceptance limits; that latter
# use is a new criterion, not a predeclared physical gate. The separate,
# tighter python_limit values are the existing Python-to-extrapolate gate.
CASES = (
    ("float", 1, (.001, .003), (.0002, .001)),
    ("spar", 3, (.0003, .0002), (.00002, .0002)),
)
if SELECTED:
    assert SELECTED in {case[0] for case in CASES}
    CASES = tuple(case for case in CASES if case[0] == SELECTED)


def _max_gap(actual, expected, limit, label):
    assert actual.shape == expected.shape, label
    gap = float(np.max(np.abs(actual - expected)))
    print(f"{label}: {gap:.9g} (limit {limit:.9g})")
    assert gap < limit, f"{label}: {gap:.9g} >= {limit:.9g}"
    return gap


def _source(name, step, checked):
    if FRESH:
        values = np.loadtxt(Path(FRESH) / name / f"{step}.csv",
                            delimiter=",", ndmin=2)
        # Fresh MATLAB must reproduce the checked source before it becomes
        # the physical target on that host.
        _max_gap(values, checked, 1e-6, f"{name} {step} source repeat")
        return values
    return checked


def _check_generated_hydro():
    if not FRESH_HYDRO:
        return
    active = (
        "properties/cg",
        "properties/disp_vol",
        "hydro_coeffs/linear_restoring_stiffness",
        "hydro_coeffs/added_mass/inf_freq",
        "hydro_coeffs/radiation_damping/impulse_response_fun/K",
        "hydro_coeffs/radiation_damping/impulse_response_fun/t",
    )
    with h5py.File(HYDRO) as bundled, h5py.File(FRESH_HYDRO) as generated:
        for body in ("body1", "body2"):
            for suffix in active:
                key = f"{body}/{suffix}"
                if suffix.endswith("impulse_response_fun/K"):
                    # Re-running BEMIO changes a few floating-point FFT
                    # round-off bits (observed max 8.19e-12), while the
                    # physical radiation kernel remains the same.
                    np.testing.assert_allclose(
                        bundled[key], generated[key], rtol=0, atol=1e-10,
                    )
                else:
                    np.testing.assert_array_equal(bundled[key], generated[key])


def _check_body_axes(case, body_name, body, source, source_trace,
                     source_limit, stride, label):
    assert source.shape == (len(source_trace), 13)
    np.testing.assert_allclose(source[:, 0], source_trace[:, 0],
                               rtol=0, atol=1e-9)
    heave_column = 1 if body_name == "float" else 3
    np.testing.assert_allclose(source[:, 3], source_trace[:, heave_column],
                               rtol=0, atol=1e-9)
    np.testing.assert_allclose(source[:, 9], source_trace[:, heave_column + 1],
                               rtol=0, atol=1e-9)
    for axis in range(6):
        if axis == 2:
            limits = source_limit if body_name == case else (1e-8, 1e-8)
        elif axis == 0:
            limits = (3e-5, 1e-5)
        elif axis == 4:
            limits = (1e-6, 1e-6)
        else:
            limits = (1e-9, 1e-9)
        _max_gap(body.position[::stride, axis], source[:, 1 + axis],
                 limits[0], f"{case} {body_name} {label} axis {axis} position")
        _max_gap(body.velocity[::stride, axis], source[:, 7 + axis],
                 limits[1], f"{case} {body_name} {label} axis {axis} speed")


@pytest.mark.parametrize("name,column,source_limit,python_limit", CASES)
def test_corrected_full_duration_source_and_independent_motion(
    name, column, source_limit, python_limit,
):
    """Advance 30 s from the published initial state without source feedback."""
    with (np.load(FIXTURES / "rm3_pto_extension_00625.npz",
                  allow_pickle=False) as old,
          np.load(FIXTURES / "rm3_pto_extension_003125.npz",
                  allow_pickle=False) as prior,
          np.load(PRIOR, allow_pickle=False) as local,
          np.load(CORRECTED, allow_pickle=False) as adopted,
          np.load(BODY_AXES, allow_pickle=False) as body_axes):
        c = old[f"{name}_0p00625"]
        d = prior[f"{name}_0p003125"]
        e = _source(name, "0p0015625", local[f"{name}_0p0015625"])
        f = _source(name, "0p00078125", local[f"{name}_0p00078125"])
        reference = adopted[name]
        assert reference.shape == (301, 9)
        # The adopted baseline retains every published 0.1 s reporting
        # sample from the exact fine-step source, without interpolation.
        np.testing.assert_allclose(
            reference, local[f"{name}_0p00078125"][::128],
            rtol=0, atol=1e-12,
        )
        for dt, values in ((.00625, c), (.003125, d),
                           (.0015625, e), (.00078125, f)):
            assert values.shape == (round(30 / dt) + 1, 9)
            np.testing.assert_allclose(values[:, 0],
                                       np.arange(len(values)) * dt,
                                       rtol=0, atol=1e-9)
            assert np.isfinite(values).all()
            assert np.max(np.abs(values[:, 7:9])) < 1e-7

        if FRESH:
            settings = np.loadtxt(Path(FRESH) / name / "settings.csv",
                                  delimiter=",", ndmin=2)
            np.testing.assert_allclose(
                settings, [[.0015625, .0015625, .0015625, .0015625,
                            .0015625, 30],
                           [.00078125, .00078125, .00078125, .00078125,
                            .00078125, 30]],
                rtol=0, atol=1e-10,
            )
        _check_generated_hydro()

        channel = slice(column, column + 2)
        old_change = np.max(np.abs(c[:, channel] - d[::2, channel]), axis=0)
        mid_change = np.max(np.abs(d[:, channel] - e[::2, channel]), axis=0)
        final_change = np.max(np.abs(e[:, channel] - f[::2, channel]), axis=0)
        assert np.all(mid_change < old_change)
        assert np.all(final_change < mid_change)
        assert np.all(final_change < source_limit), (name, final_change)

        # Preserve the earlier convergence diagnostic at every common
        # 0.00625 s sample. The direct corrected-source gate is below.
        old_extrap = 2 * d[::2, channel] - c[:, channel]
        mid_extrap = 2 * e[::4, channel] - d[::2, channel]
        final_extrap = 2 * f[::8, channel] - e[::4, channel]
        old_spread = np.max(np.abs(mid_extrap - old_extrap), axis=0)
        final_spread = np.max(np.abs(final_extrap - mid_extrap), axis=0)
        assert np.all(final_spread < old_spread), (name, final_spread)

        wec = WEC(f"RM3 PTO extension {name} corrected validation")
        float_body = wec.body("float", HYDRO,
                              inertia=(0, 21_306_090.66, 0))
        spar = wec.body("spar", HYDRO,
                        inertia=(0, 94_407_091.24, 0))
        wec.floating_joint(float_body, spar, added_mass_scheme="implicit")
        initial = {"float_heave": 5} if name == "float" else {"spar_heave": -5}
        result = wec.run(NoWave(), dt=.00078125, end_time=30,
                         radiation_memory=60, initial_coordinate=initial)
        assert result.wave_elevation is None
        np.testing.assert_allclose(result.time, f[:, 0], rtol=0, atol=1e-9)

        active = result.bodies[name]
        inactive = result.bodies["spar" if name == "float" else "float"]
        _max_gap(active.position[:, 2], f[:, column], source_limit[0],
                 f"{name} active body heave")
        _max_gap(active.velocity[:, 2], f[:, column + 1], source_limit[1],
                 f"{name} active body heave speed")
        inactive_column = 3 if name == "float" else 1
        _max_gap(inactive.position[:, 2], f[:, inactive_column], 1e-8,
                 f"{name} inactive body heave")
        _max_gap(inactive.velocity[:, 2], f[:, inactive_column + 1], 1e-8,
                 f"{name} inactive body heave speed")
        # Six-axis comparison is portable at all 301 adopted reporting
        # samples. Fresh MATLAB additionally checks every fine-grid sample;
        # small coupled surge is compared source-relative, not to zero.
        for body_name in ("float", "spar"):
            body = result.bodies[body_name]
            assert np.isfinite(body.position).all()
            assert np.isfinite(body.velocity).all()
            _check_body_axes(
                name, body_name, body, body_axes[f"{name}_{body_name}"],
                reference, source_limit, 128, "adopted",
            )
            if FRESH:
                source_body = np.loadtxt(
                    Path(FRESH) / name / f"0p00078125_{body_name}_body.csv",
                    delimiter=",", ndmin=2,
                )
                _check_body_axes(
                    name, body_name, body, source_body, f,
                    source_limit, 1, "fresh",
                )

        pto = result.ptos["relative_heave"]
        _max_gap(pto.stroke, f[:, 5], source_limit[0], f"{name} PTO stroke")
        _max_gap(pto.velocity, f[:, 6], source_limit[1], f"{name} PTO speed")
        _max_gap(pto.force, f[:, 7], 1e-7, f"{name} zero PTO force")
        _max_gap(pto.force * pto.velocity, f[:, 8], 1e-7,
                 f"{name} zero PTO power")

        # Check the independently advanced public result against the
        # repository-owned 301-sample corrected validation baseline too.
        sampled = np.column_stack((
            result.time[::128],
            result.bodies["float"].position[::128, 2],
            result.bodies["float"].velocity[::128, 2],
            result.bodies["spar"].position[::128, 2],
            result.bodies["spar"].velocity[::128, 2],
            pto.stroke[::128], pto.velocity[::128],
            pto.force[::128],
            (pto.force * pto.velocity)[::128],
        ))
        for index, limit in ((1, source_limit[0]), (2, source_limit[1]),
                             (3, source_limit[0]), (4, source_limit[1]),
                             (5, source_limit[0]), (6, source_limit[1]),
                             (7, 1e-7), (8, 1e-7)):
            if index in ((3, 4) if name == "float" else (1, 2)):
                limit = 1e-8
            _max_gap(sampled[:, index], reference[:, index], limit,
                     f"{name} adopted baseline column {index}")

        python = np.column_stack((active.position[::8, 2],
                                  active.velocity[::8, 2]))
        extrapolation_gap = np.max(np.abs(python - final_extrap), axis=0)
        assert np.all(extrapolation_gap < final_spread), (
            name, extrapolation_gap, final_spread,
        )
        assert np.all(extrapolation_gap < python_limit), (
            name, extrapolation_gap, python_limit,
        )

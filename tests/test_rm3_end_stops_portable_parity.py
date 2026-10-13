"""Run the full refined RM3 End_Stops case against pinned MATLAB in CI."""

from pathlib import Path

import numpy as np

from wecsim import LinearHardStops, solve_rm3_regular


ROOT = Path(__file__).resolve().parents[1]
HYDRO = ROOT / "examples/data/rm3_end_stops.h5"
MATLAB = ROOT / "tests/fixtures/matlab/rm3_end_stops_full.npz"
AXES = (0, 2, 4)


def _bound(actual, expected, limit, label):
    assert actual.shape == expected.shape, label
    error = float(np.max(np.abs(actual - expected)))
    assert np.isfinite(error) and error < limit, (
        f"{label}: {error:.8g} >= {limit:.8g}"
    )


def test_full_duration_refined_matlab_end_stops_are_stable():
    """Keep the source-step check alongside the independent Python gate."""
    with np.load(MATLAB, allow_pickle=False) as source:
        for label, dt, count in (("dt0025", .025, 16001),
                                 ("dt00125", .0125, 32001)):
            np.testing.assert_allclose(source[f"{label}_run_settings"],
                                       (dt, 400, 100), rtol=0, atol=1e-12)
            np.testing.assert_allclose(source[f"{label}_settings"],
                                       (-.6, .6, 1e8, 1e8, 0, 0, .0001, .0001),
                                       rtol=0, atol=1e-12)
            np.testing.assert_allclose(source[f"{label}_time"],
                                       np.arange(count) * dt, rtol=0,
                                       atol=1e-10)
            assert source[f"{label}_body_position"].shape == (count, 2, 3)
            assert source[f"{label}_body_velocity"].shape == (count, 2, 3)
            assert source[f"{label}_pto"].shape == (count, 3)

        for kind, limits in (
            ("body_position", (.0025, .004, 8e-5)),
            ("body_velocity", (.004, .022, 1e-4)),
        ):
            coarse = source[f"dt0025_{kind}"]
            fine = source[f"dt00125_{kind}"][::2]
            for body in range(2):
                for axis, limit in enumerate(limits):
                    _bound(coarse[:, body, axis], fine[:, body, axis],
                           limit, f"source body {body + 1} {kind} {AXES[axis]}")
        coarse, fine = source["dt0025_pto"], source["dt00125_pto"][::2]
        for axis, limit in enumerate((.004, .025, 350_000)):
            _bound(coarse[:, axis], fine[:, axis], limit,
                   f"source PTO channel {axis}")
        coarse_energy = np.trapezoid(
            1_200_000 * coarse[:, 1]**2, source["dt0025_time"],
        )
        fine_energy = np.trapezoid(
            1_200_000 * source["dt00125_pto"][:, 1]**2,
            source["dt00125_time"],
        )
        assert fine_energy > 0
        assert abs(coarse_energy - fine_energy) / fine_energy < .006


def test_independent_full_duration_rm3_end_stops_motion():
    """Advance 400 s from HDF5 and stop law, never source motion or force."""
    response = solve_rm3_regular(
        HYDRO, pto_hard_stops=LinearHardStops(-.6, .6, 1e8, 1e8),
        dt=.0125, end_time=400,
    )
    with np.load(MATLAB, allow_pickle=False) as source:
        np.testing.assert_allclose(response.time, source["dt00125_time"],
                                   rtol=0, atol=1e-10)
        for kind, limits in (
            ("body_position", (.005, .006, 1e-4)),
            ("body_velocity", (.005, .03, 1e-4)),
        ):
            actual = getattr(response, kind)[:, :, AXES]
            expected = source[f"dt00125_{kind}"]
            for body in range(2):
                for axis, limit in enumerate(limits):
                    _bound(actual[:, body, axis], expected[:, body, axis],
                           limit, f"Python body {body + 1} {kind} {AXES[axis]}")
        pto = source["dt00125_pto"]
        for actual, axis, limit, label in (
            (response.pto_stroke, 0, .006, "PTO stroke"),
            (response.pto_velocity, 1, .035, "PTO speed"),
            (response.pto_force, 2, 450_000, "PTO force"),
        ):
            _bound(actual, pto[:, axis], limit, label)
        python_energy = np.trapezoid(
            response.pto_dissipated_power, response.time,
        )
        matlab_energy = np.trapezoid(1_200_000 * pto[:, 1]**2,
                                     source["dt00125_time"])
        assert abs(python_energy - matlab_energy) / matlab_energy < .007

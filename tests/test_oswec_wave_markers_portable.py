"""Keep the published OSWEC marker and motion pair active in Python CI."""

from pathlib import Path

import numpy as np

from wecsim import RegularWave, WEC, WorldPoint


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "tests/fixtures/matlab/oswec_wave_markers.npz"
HYDRO = ROOT / "examples/data/oswec.h5"


def _max_error(actual, expected, limit, label):
    assert actual.shape == expected.shape, f"{label}: mismatched shape"
    assert np.isfinite(actual).all() and np.isfinite(expected).all(), label
    error = float(np.max(np.abs(actual - expected)))
    assert error < limit, f"{label}: {error:.6g} exceeds {limit:.6g}"
    return error


def test_published_oswec_marker_elevations():
    """Pair all 25 markers and origin elevation at every published time."""
    with np.load(SOURCE, allow_pickle=False) as source:
        time = source["time"]
        locations = source["locations"]
        height, period, direction, ramp, gravity, depth, source_k = (
            source["wave_parameters"]
        )
        np.testing.assert_array_equal(
            source["wave_parameters"][:6], [2.5, 8, 0, 100, 9.81, 10.9],
        )
        assert time.shape == (4001,)
        np.testing.assert_allclose(time, np.arange(4001) * .1,
                                   rtol=0, atol=1e-10)
        grid = np.arange(-20, 21, 10)
        x, y = np.meshgrid(grid, grid)
        np.testing.assert_array_equal(
            locations, np.column_stack((x.ravel(order="F"),
                                        y.ravel(order="F"))),
        )
        actual = RegularWave(height, period, direction).elevation_at(
            time, locations, water_depth=depth, ramp_time=ramp, g=gravity,
        )
        _max_error(actual, source["markers"], 1e-10, "OSWEC markers")
        center = np.flatnonzero(np.all(locations == [0, 0], axis=1))
        assert center.size == 1
        _max_error(actual[:, center[0]], source["origin"], 1e-10,
                   "OSWEC origin elevation")
        omega = 2 * np.pi / period
        dispersion_error = abs(
            omega**2 - gravity * source_k * np.tanh(source_k * depth)
        ) / omega**2
        assert dispersion_error < 1e-10


def test_published_oswec_wave_markers_motion_and_pto():
    """Gate both body histories and all hinge PTO channels for 400 s."""
    wec = WEC("OSWEC published Wave Markers")
    flap = wec.body("flap", HYDRO, mass=127_000, inertia=(0, 1.85e6, 0))
    base = wec.body("base", HYDRO, mass=999, inertia=(999, 999, 999))
    wec.fixed_hinge(flap, base, location=WorldPoint(0, 0, -10),
                    pto_location=WorldPoint(0, 0, -8.9), damping=12_000)
    result = wec.run(RegularWave(2.5, 8), dt=.1, end_time=400, ramp_time=100)

    with np.load(SOURCE, allow_pickle=False) as source:
        np.testing.assert_allclose(result.time, source["time"],
                                   rtol=0, atol=1e-10)
        assert source["body_position"].shape == (2, 4001, 6)
        assert source["body_velocity"].shape == (2, 4001, 6)
        position_limits = {0: .03, 2: .012, 4: .007}
        speed_limits = {0: .025, 2: .01, 4: .006}
        for index, name in enumerate(("flap", "base")):
            for dof in range(6):
                position_limit = (position_limits.get(dof, 1e-10)
                                  if index == 0 else 1e-10)
                speed_limit = (speed_limits.get(dof, 1e-10)
                               if index == 0 else 1e-10)
                _max_error(result.bodies[name].position[:, dof],
                           source["body_position"][index, :, dof],
                           position_limit, f"{name} position DOF {dof}")
                _max_error(result.bodies[name].velocity[:, dof],
                           source["body_velocity"][index, :, dof],
                           speed_limit, f"{name} velocity DOF {dof}")
        pto = result.ptos["hinge"]
        for label, actual, expected, limit in (
            ("stroke", pto.stroke, source["pto_stroke"], .006),
            ("speed", pto.velocity, source["pto_speed"], .005),
            ("force", pto.force, source["pto_force"], 60),
            # MATLAB reports mechanical power with the opposite sign.
            ("absorbed power", pto.absorbed_power, -source["pto_power"], 30),
        ):
            _max_error(actual, expected, limit, f"OSWEC PTO {label}")

        # The passive damper's work uses velocity squared; the separate
        # source-step audit checks cumulative energy over the full case.
        np.testing.assert_allclose(
            pto.absorbed_power, 12_000 * pto.velocity**2,
            rtol=0, atol=1e-8,
        )

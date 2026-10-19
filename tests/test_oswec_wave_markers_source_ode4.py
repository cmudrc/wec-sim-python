"""Gate the public opt-in solver against published OSWEC regular-wave motion."""

from pathlib import Path

import numpy as np

from wecsim import RegularWave, WEC, WorldPoint


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "tests/fixtures/matlab/oswec_wave_markers.npz"
# The published application generates this 36-heading BEMIO file. The bundled
# example HDF5 has only 10 headings and slightly different coefficients.
HYDRO = ROOT / "tests/fixtures/matlab/oswec_wave_markers/oswec.h5"


def _max_error(actual, expected, limit, label):
    assert actual.shape == expected.shape, label
    assert np.isfinite(actual).all() and np.isfinite(expected).all(), label
    error = float(np.max(np.abs(actual - expected)))
    assert error < limit, f"{label}: {error:.9g} exceeds {limit:.9g}"


def test_full_duration_public_source_ode4_motion_pto_and_work():
    """Advance from HDF5 and wave inputs, without replaying source states."""
    wec = WEC("OSWEC Wave Markers published-step source solver")
    flap = wec.body("flap", HYDRO, mass=127_000, inertia=(0, 1.85e6, 0))
    base = wec.body("base", HYDRO, mass=999, inertia=(999, 999, 999))
    wec.fixed_hinge(
        flap, base, location=WorldPoint(0, 0, -10),
        pto_location=WorldPoint(0, 0, -8.9), damping=12_000,
        added_mass_scheme="simulink_delay",
    )
    result = wec.run(
        RegularWave(2.5, 8), dt=.1, end_time=400, ramp_time=100,
        solver="wecsim_ode4",
    )

    with np.load(SOURCE, allow_pickle=False) as source:
        time = source["time"]
        np.testing.assert_allclose(result.time, time, rtol=0, atol=1e-10)
        assert time.shape == (4001,)
        for index, name in enumerate(("flap", "base")):
            for axis in range(6):
                position_limit = ({0: 1e-5, 2: 2e-6, 4: 5e-6}.get(axis, 1e-10)
                                  if index == 0 else 1e-10)
                speed_limit = ({0: 5e-6, 2: 2e-6, 4: 5e-6}.get(axis, 1e-10)
                               if index == 0 else 1e-10)
                _max_error(result.bodies[name].position[:, axis],
                           source["body_position"][index, :, axis],
                           position_limit, f"{name} position DOF {axis}")
                _max_error(result.bodies[name].velocity[:, axis],
                           source["body_velocity"][index, :, axis],
                           speed_limit, f"{name} velocity DOF {axis}")

        pto = result.ptos["hinge"]
        for label, actual, expected, limit in (
            ("stroke", pto.stroke, source["pto_stroke"], 5e-6),
            ("speed", pto.velocity, source["pto_speed"], 5e-6),
            ("torque", pto.force, source["pto_force"], .1),
            ("absorbed power", pto.absorbed_power, -source["pto_power"], .1),
        ):
            _max_error(actual, expected, limit, f"PTO {label}")

        for label, start, stop in (
            ("ramp", 0, 100),
            ("early steady", 100, 250),
            ("late steady", 250, 400),
            ("full", 0, 400),
        ):
            window = (time >= start) & (time <= stop)
            source_work = -np.trapezoid(source["pto_power"][window], time[window])
            python_work = np.trapezoid(pto.absorbed_power[window], time[window])
            assert source_work > 0, label
            relative_error = abs(python_work - source_work) / source_work
            assert relative_error < 1e-4, (label, relative_error)

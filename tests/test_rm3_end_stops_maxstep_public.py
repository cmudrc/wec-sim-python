"""Gate public RM3 End_Stops motion against a MaxStep-only MATLAB arm."""

from pathlib import Path

import numpy as np

from wecsim import LinearHardStops, RegularWave, WEC


ROOT = Path(__file__).resolve().parents[1]
HYDRO = ROOT / "examples/data/rm3_end_stops.h5"
SOURCE = ROOT / "tests/fixtures/matlab/rm3_end_stops_maxstep_120s.npz"


def test_independent_public_motion_against_maxstep_only_source():
    """Advance the physical Python default from HDF5, without source states."""
    with np.load(SOURCE, allow_pickle=False) as source:
        np.testing.assert_allclose(source["settings"], (.1, 120, 100, .025),
                                   rtol=0, atol=1e-12)
        time = source["time"]
        wave = source["wave"]
        source_pto = source["pto"]
        source_position = source["body_position"]
        source_velocity = source["body_velocity"]
    np.testing.assert_allclose(time, np.arange(1201) * .1,
                               rtol=0, atol=1e-10)
    np.testing.assert_allclose(wave[:, 0], time, rtol=0, atol=1e-10)
    assert source_pto.shape == (1201, 3)
    assert source_position.shape == source_velocity.shape == (1201, 2, 3)
    assert all(np.isfinite(array).all() for array in (
        wave, source_pto, source_position, source_velocity,
    ))

    wec = WEC("RM3 End Stops MaxStep audit")
    float_body = wec.body("float", HYDRO,
                          inertia=(0, 21_306_090.66, 0))
    spar = wec.body("spar", HYDRO,
                    inertia=(0, 94_407_091.24, 0))
    wec.floating_joint(
        float_body, spar, damping=1_200_000,
        hard_stops=LinearHardStops(-.6, .6, 1e8, 1e8),
    )
    result = wec.run(RegularWave(2.5, 8), dt=.1, end_time=120,
                     ramp_time=100)
    np.testing.assert_allclose(result.time, time, rtol=0, atol=1e-9)
    for index, name in enumerate(("float", "spar")):
        body = result.bodies[name]
        for kind, actual, expected, limits in (
            ("position", body.position[:, (0, 2, 4)],
             source_position[:, index], (.005, .006, 1e-4)),
            ("velocity", body.velocity[:, (0, 2, 4)],
             source_velocity[:, index], (.005, .03, 1e-4)),
        ):
            for axis, limit in enumerate(limits):
                error = np.max(np.abs(actual[:, axis] - expected[:, axis]))
                print(f"{name} {kind} axis {axis}: {error:.9g}")
                assert error < limit, f"{name} {kind} {axis}: {error} >= {limit}"

    pto = result.ptos["relative_heave"]
    for actual, expected, limit, name in (
        (pto.stroke, source_pto[:, 0], .006, "PTO stroke"),
        (pto.velocity, source_pto[:, 1], .035, "PTO speed"),
        (pto.force, source_pto[:, 2], 450_000, "PTO force"),
    ):
        error = np.max(np.abs(actual - expected))
        print(f"{name}: {error:.9g}")
        assert error < limit, f"{name}: {error} >= {limit}"
    python_work = np.trapezoid(pto.absorbed_power, result.time)
    source_work = np.trapezoid(1_200_000 * source_pto[:, 1] ** 2, time)
    assert source_work > 0
    work_error = abs(python_work - source_work) / source_work
    print(f"relative damper work: {work_error:.9g}")
    assert work_error < .007

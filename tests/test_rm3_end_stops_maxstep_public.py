"""Gate public RM3 End_Stops motion against a MaxStep-only MATLAB arm."""

from pathlib import Path

import numpy as np
import pytest

from wecsim import LinearHardStops, RegularWave, WEC


ROOT = Path(__file__).resolve().parents[1]
HYDRO = ROOT / "examples/data/rm3_end_stops.h5"
SOURCE_120 = ROOT / "tests/fixtures/matlab/rm3_end_stops_maxstep_120s.npz"
SOURCE_400 = ROOT / "tests/fixtures/matlab/rm3_end_stops_maxstep_400s.npz"
REFINED = ROOT / "tests/fixtures/matlab/rm3_end_stops_full.npz"


@pytest.mark.parametrize("source_path", (SOURCE_120, SOURCE_400),
                         ids=("120s", "400s"))
def test_independent_public_motion_against_maxstep_only_source(source_path):
    """Advance the physical Python default from HDF5, without source states."""
    with np.load(source_path, allow_pickle=False) as source:
        end_time = int(source["settings"][1])
        assert end_time in (120, 400)
        np.testing.assert_allclose(source["settings"],
                                   (.1, end_time, 100, .025),
                                   rtol=0, atol=1e-12)
        time = source["time"]
        wave = source["wave"]
        source_pto = source["pto"]
        source_position = source["body_position"]
        source_velocity = source["body_velocity"]
    count = round(end_time / .1) + 1
    np.testing.assert_allclose(time, np.arange(count) * .1,
                               rtol=0, atol=1e-10)
    np.testing.assert_allclose(wave[:, 0], time, rtol=0, atol=1e-10)
    assert source_pto.shape == (count, 3)
    assert source_position.shape == source_velocity.shape == (count, 2, 3)
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
    result = wec.run(RegularWave(2.5, 8), dt=.1, end_time=end_time,
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


def test_full_duration_maxstep_source_closes_published_step_gap():
    """Preserve the MATLAB self-change evidence behind the derived gate."""
    with (np.load(SOURCE_120, allow_pickle=False) as prefix,
          np.load(SOURCE_400, allow_pickle=False) as source,
          np.load(REFINED, allow_pickle=False) as refined):
        np.testing.assert_allclose(source["pto"][:1201], prefix["pto"],
                                   rtol=0, atol=1e-2)
        np.testing.assert_allclose(
            source["body_position"][:1201], prefix["body_position"],
            rtol=0, atol=1e-7,
        )
        np.testing.assert_allclose(
            source["body_velocity"][:1201], prefix["body_velocity"],
            rtol=0, atol=1e-6,
        )
        time = source["time"]
        fine = refined["dt00125_pto"][::8]
        assert source["published_pto"].shape == source["pto"].shape == fine.shape
        after_contact = time >= 68.0
        errors = []
        for pto in (source["published_pto"], source["pto"]):
            stroke_gap = pto[:, 0] - fine[:, 0]
            rms = np.sqrt(np.mean(stroke_gap[after_contact] ** 2))
            maximum = np.max(np.abs(stroke_gap[after_contact]))
            work = np.trapezoid(1_200_000 * pto[:, 1] ** 2, time)
            fine_work = np.trapezoid(1_200_000 * fine[:, 1] ** 2, time)
            errors.append((rms, maximum, abs(work - fine_work) / fine_work))
    for improved, original, fraction in zip(
            errors[1], errors[0], (.1, .1, .25), strict=True):
        assert original > 0
        assert improved < fraction * original

"""Full-duration public RM3 MooringMatrix gate using pinned MATLAB inputs."""

from pathlib import Path

import numpy as np

from wecsim import ImportedElevationWave, WEC


FIXTURE = Path(__file__).parent / "fixtures/matlab/rm3_mooring_matrix"


def _within(actual, expected, limit, label):
    actual = np.asarray(actual)
    expected = np.asarray(expected)
    assert actual.shape == expected.shape, label
    assert np.isfinite(actual).all() and np.isfinite(expected).all(), label
    error = float(np.max(np.abs(actual - expected)))
    assert error < limit, f"{label}: {error:.6g} exceeds {limit:.6g}"


def test_published_rm3_imported_sea_joint_and_pto_over_400_seconds():
    with np.load(FIXTURE / "source.npz") as archive:
        source = {name: archive[name] for name in archive.files}
    assert source["time"].shape == (40001,)
    np.testing.assert_allclose(source["time"], np.arange(40001) * .01,
                               rtol=0, atol=1e-9)

    wec = WEC("RM3 MooringMatrix")
    float_body = wec.body("float", FIXTURE / "rm3.h5",
                          inertia=(0, 21_306_090.66, 0))
    spar = wec.body("spar", FIXTURE / "rm3.h5",
                    inertia=(0, 94_407_091.24, 0))
    wec.floating_joint(float_body, spar, damping=1_200_000,
                       mooring_surge_stiffness=100_000)
    result = wec.run(
        ImportedElevationWave(FIXTURE / "etaData.mat",
                              reapply_force_ramp=True),
        dt=.01, end_time=400, ramp_time=40, radiation_memory=60,
        initial_coordinate={"spar_heave": -.21},
    )
    np.testing.assert_allclose(result.time, source["time"], rtol=0, atol=1e-9)
    _within(result.wave_elevation, source["wave"], 1e-10, "imported sea")

    for number, name in enumerate(("float", "spar")):
        body = result.bodies[name]
        for column, label, position_limit, speed_limit in (
            (0, "surge", .015, .006),
            (2, "heave", .004, .0006),
            (4, "pitch", .0009, .0004),
        ):
            _within(body.position[:, column],
                    source["body_active_position"][:, number, column // 2],
                    position_limit, f"{name} {label} position")
            _within(body.velocity[:, column],
                    source["body_active_velocity"][:, number, column // 2],
                    speed_limit, f"{name} {label} speed")

    pto = result.ptos["relative_heave"]
    _within(pto.stroke - pto.stroke[0], source["pto_stroke"],
            .004, "PTO stroke from initial datum")
    _within(pto.velocity, source["pto_speed"], .0005, "PTO speed")
    _within(pto.force, source["pto_force"], 500, "PTO force")
    _within(pto.absorbed_power, -source["pto_power"], 500,
            "PTO absorbed power")
    source_work = -np.trapezoid(source["pto_power"], x=source["time"])
    python_work = np.trapezoid(pto.absorbed_power, x=result.time)
    relative_work_gap = abs(python_work - source_work) / source_work
    print(f"RM3 MooringMatrix relative PTO work gap: {relative_work_gap:.6%}")
    assert relative_work_gap < .01

    mooring = dict(result.raw.extra_outputs)
    _within(mooring["mooring_surge_position"],
            source["mooring_surge_position"], .015,
            "mooring surge position")
    _within(mooring["mooring_surge_force"],
            source["mooring_surge_force"], 1_500,
            "mooring surge force")

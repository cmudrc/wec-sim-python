"""Run the published 400 s directional OSWEC case in ordinary Python CI."""

from hashlib import sha256
from pathlib import Path

import numpy as np

from wecsim import FullDirectionalSpectrumWave, WEC, WorldPoint
from wecsim.environment.irregularWave import (
    imported_full_directional_components, synthesize_full_directional_response,
)


FIXTURE = Path(__file__).parent / "fixtures/matlab/oswec_full_directional"
HYDRO = Path(__file__).parent / "fixtures/matlab/oswec_multi_wave/oswec.h5"
SPECTRUM = FIXTURE / "fullDirSpectrum.mat"
SOURCE = FIXTURE / "source.npz"
HYDRO_SHA256 = "b9ccf3a32e0893c23b1645c15a2baef4a0c6f678b90d177d5797e3c65de4b7d3"
SPECTRUM_SHA256 = "c0abee4fe5da352412d64342ea276832bfb19595774a0f1f2106f5f109dfaee9"
SOURCE_SHA256 = "52a83e9d86b18b16ea876d6d8ea88cd6c04db80b189ba1717d3b661fac6e121b"


def _within(actual, expected, limit, label):
    actual = np.asarray(actual)
    expected = np.asarray(expected)
    assert actual.shape == expected.shape, label
    assert np.isfinite(actual).all() and np.isfinite(expected).all(), label
    error = float(np.max(np.abs(actual - expected)))
    assert error < limit, f"{label}: {error:.6g} exceeds {limit:.6g}"


def test_published_full_directional_wave_forces_motion_and_pto():
    for path, digest in ((HYDRO, HYDRO_SHA256),
                         (SPECTRUM, SPECTRUM_SHA256),
                         (SOURCE, SOURCE_SHA256)):
        assert sha256(path.read_bytes()).hexdigest() == digest
    with np.load(SOURCE) as archive:
        source = {name: archive[name] for name in archive.files}
    assert source["phase"].shape == (47, 180)
    assert source["body_excitation"].shape == (8001, 2, 6)
    np.testing.assert_allclose(source["time"], np.arange(8001) * .05,
                               rtol=0, atol=1e-10)

    components = imported_full_directional_components(
        HYDRO, SPECTRUM, seed=1, phase_generator="matlab",
    )
    _within(components.phase, source["phase"], 1e-13, "MATLAB seed phases")
    _within(components.d_theta, np.full(180, np.deg2rad(2)),
            1e-12, "heading width")
    for number in (1, 2):
        incident = synthesize_full_directional_response(
            HYDRO, components, dt=.05, end_time=400, ramp_time=100,
            body_number=number, excitation_interpolation="spline_frequency",
            force_quadrature="matlab_omitted",
        )
        _within(incident.elevation, source["elevation"], 1e-11,
                "published wave elevation")
        _within(incident.excitation_force,
                source["body_excitation"][:, number - 1], 1e-6,
                f"body {number} published excitation")

    # The source force block omits the heading-bin width. Keep the physical
    # quadrature as the default while explicitly reproducing that source path.
    physical = synthesize_full_directional_response(
        HYDRO, components, dt=.05, end_time=400, ramp_time=100,
        body_number=1, excitation_interpolation="spline_frequency",
    )
    _within(physical.excitation_force,
            source["body_excitation"][:, 0] * np.sqrt(np.deg2rad(2)),
            1e-6, "physical heading quadrature")

    wec = WEC("Published full directional OSWEC")
    flap = wec.body("flap", HYDRO, mass=127000,
                    inertia=(0, 1.85e6, 0), hydro_body=1)
    base = wec.body("base", HYDRO, mass=999,
                    inertia=(999, 999, 999), hydro_body=2)
    wec.fixed_hinge(
        flap, base, location=WorldPoint(0, 0, -10),
        pto_location=WorldPoint(0, 0, -8.9), damping=12000,
        added_mass_scheme="simulink_delay",
    )
    result = wec.run(
        FullDirectionalSpectrumWave(
            SPECTRUM, seed=1, phase_generator="matlab",
            excitation_interpolation="spline_frequency",
            force_quadrature="matlab_omitted",
        ),
        dt=.05, end_time=400, ramp_time=100, radiation_memory=30,
    )
    np.testing.assert_allclose(result.time, source["time"], rtol=0, atol=1e-10)
    _within(result.wave_elevation, source["elevation"], 1e-11,
            "public wave elevation")

    for number, name in enumerate(("flap", "base")):
        body = result.bodies[name]
        if name == "base":
            np.testing.assert_allclose(
                body.position, source["body_position"][:, number],
                rtol=0, atol=1e-10,
            )
            np.testing.assert_allclose(
                body.velocity, source["body_velocity"][:, number],
                rtol=0, atol=1e-10,
            )
            continue
        for dof, position_limit, speed_limit in (
            (0, .002, .003), (2, .0004, .0007), (4, .0004, .0006),
        ):
            _within(body.position[:, dof],
                    source["body_position"][:, number, dof],
                    position_limit, f"{name} position {dof}")
            _within(body.velocity[:, dof],
                    source["body_velocity"][:, number, dof],
                    speed_limit, f"{name} speed {dof}")
        for dof in (1, 3, 5):
            np.testing.assert_allclose(
                body.position[:, dof], source["body_position"][:, number, dof],
                rtol=0, atol=1e-10,
            )
            np.testing.assert_allclose(
                body.velocity[:, dof], source["body_velocity"][:, number, dof],
                rtol=0, atol=1e-10,
            )

    raw = result.raw
    _within(raw.pto_force, source["pto_force"], 8,
            "published PTO torque")
    power = raw.pto_force * raw.body_velocity[:, 0, 4]
    _within(power, source["pto_power"], 5,
            "published PTO power")
    source_work = np.trapezoid(source["pto_power"], x=source["time"])
    python_work = np.trapezoid(power, x=result.time)
    assert source_work < 0
    assert abs(python_work - source_work) / abs(source_work) < .01

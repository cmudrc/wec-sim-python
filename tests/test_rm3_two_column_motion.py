"""Pair a seeded two-column imported sea with independent RM3 motion."""

from hashlib import sha256
from pathlib import Path

import numpy as np
from scipy.integrate import trapezoid

from wecsim import ImportedSpectrumWave, WEC
from wecsim.environment.irregularWave import imported_spectrum_components


ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests/fixtures/matlab/rm3_two_column_motion/source.npz"
SPECTRUM = ROOT / "tests/fixtures/matlab/imported_spectrum_two_column/spectrum.mat"
HYDRO = ROOT / "examples/data/rm3.h5"


def _max_error(actual, expected, limit, label):
    assert actual.shape == expected.shape, label
    assert np.isfinite(actual).all() and np.isfinite(expected).all(), label
    error = float(np.max(np.abs(actual - expected)))
    assert error < limit, f"{label}: {error:.9g} >= {limit:.9g}"


def test_seeded_two_column_source_is_pinned():
    assert sha256(FIXTURE.read_bytes()).hexdigest() == (
        "4b590fa83fde281abe2ce839409c730aa6bd58b81358fe52cf5477b6e34855e9"
    )


def test_public_rm3_seeded_two_column_motion_matches_local_matlab():
    """Use only the spectrum and bundled HDF5 as Python simulation inputs."""
    with np.load(FIXTURE, allow_pickle=False) as saved:
        wave = saved["wave"]
        body1 = saved["body1"]
        body2 = saved["body2"]
        pto = saved["pto"]
        components = saved["components"]
    assert wave.shape == (301, 2)
    assert body1.shape == body2.shape == (301, 13)
    assert pto.shape == (301, 25)
    assert components.shape == (5, 2)
    np.testing.assert_allclose(body1[:, 0], wave[:, 0], rtol=0, atol=1e-10)
    np.testing.assert_allclose(body2[:, 0], wave[:, 0], rtol=0, atol=1e-10)
    np.testing.assert_allclose(pto[:, 0], wave[:, 0], rtol=0, atol=1e-10)

    sea = imported_spectrum_components(
        HYDRO, SPECTRUM, seed=1, phase_generator="matlab",
    )
    np.testing.assert_allclose(sea.omega, components[:, 0], rtol=0, atol=1e-13)
    np.testing.assert_allclose(sea.phase[:, 0], components[:, 1],
                               rtol=0, atol=1e-13)

    wec = WEC("seeded two-column RM3")
    float_body = wec.body("float", HYDRO, inertia=(0, 21_306_090.66, 0))
    spar = wec.body("spar", HYDRO, inertia=(0, 94_407_091.24, 0))
    # The source's delayed added-mass feedback is explicit for numerical
    # comparison; the public physical default remains implicit.
    wec.floating_joint(float_body, spar, damping=1_200_000,
                       added_mass_scheme="simulink_delay")
    result = wec.run(
        ImportedSpectrumWave(SPECTRUM, seed=1, phase_generator="matlab"),
        dt=.1, end_time=30, ramp_time=1, radiation_memory=60,
    )
    assert result.case["wave"]["seed"] == 1
    assert result.case["simulation"]["added_mass_scheme"] == "simulink_delay"
    np.testing.assert_allclose(result.time, wave[:, 0], rtol=0, atol=1e-10)
    _max_error(result.wave_elevation, wave[:, 1], 1e-12, "wave elevation")

    limits = {
        "surge": (.0005, .0012),
        "heave": (.0002, .0004),
        "pitch": (1e-5, 2e-5),
    }
    for name, source in (("float", body1), ("spar", body2)):
        actual = result.bodies[name]
        for label, axis in (("surge", 0), ("heave", 2), ("pitch", 4)):
            position_limit, speed_limit = limits[label]
            _max_error(actual.position[:, axis], source[:, 1 + axis],
                       position_limit, f"{name} {label} position")
            _max_error(actual.velocity[:, axis], source[:, 7 + axis],
                       speed_limit, f"{name} {label} speed")

    actual_pto = result.ptos["relative_heave"]
    for label, actual, expected, limit in (
        ("stroke", actual_pto.stroke, pto[:, 3], .00025),
        ("speed", actual_pto.velocity, pto[:, 9], .00045),
        ("force", actual_pto.force, pto[:, 15], 500),
        ("power", actual_pto.absorbed_power, -pto[:, 21], 250),
    ):
        _max_error(actual, expected, limit, f"PTO {label}")
    source_work = -trapezoid(pto[:, 21], x=pto[:, 0])
    python_work = trapezoid(actual_pto.absorbed_power, x=result.time)
    assert source_work > 0
    assert abs(python_work - source_work) / source_work < .002

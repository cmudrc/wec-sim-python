"""Keep the published floating OWC wave and excitation pair in ordinary CI."""

from hashlib import sha256
from pathlib import Path

import numpy as np

from wecsim.dynamics.floatingOwcDynamics import _body_data


FIXTURE = Path(__file__).parent / "fixtures/matlab/floating_owc_wave_force"
HYDRO = FIXTURE / "floatingOWC.h5"
SOURCE = FIXTURE / "source.npz"
HYDRO_SHA256 = "c9092b10cf35bc438805f719da5434478a2211ca9261e9a33d562dbf4d4e1a11"


def test_published_floating_owc_wave_and_both_body_excitations():
    """Rebuild 500 s of source forcing from the exact BEM file, not source poses."""
    assert sha256(HYDRO.read_bytes()).hexdigest() == HYDRO_SHA256
    with np.load(SOURCE, allow_pickle=False) as archived:
        time = archived["time"].copy()
        source_wave = archived["wave"].copy()
        source_forces = (archived["body1_excitation"].copy(),
                         archived["body2_excitation"].copy())

    np.testing.assert_allclose(time, np.arange(50_001) * .01,
                               rtol=0, atol=1e-10)
    assert source_wave.shape == (50_001,)
    assert all(force.shape == (50_001, 6) for force in source_forces)
    assert np.isfinite(source_wave).all()
    assert all(np.isfinite(force).all() for force in source_forces)

    frequency = 2 * np.pi / 11.2
    ramp = np.ones_like(time)
    early = time < 50
    ramp[early] = (1 - np.cos(np.pi * time[early] / 50)) / 2
    cosine = np.cos(frequency * time)
    sine = np.sin(frequency * time)
    wave = 2.25 * ramp * cosine
    assert np.max(np.abs(wave - source_wave)) < 1e-11

    column_mass = 4_493_450.0
    radius = 5.89 / 2
    transverse = column_mass * (3 * radius**2 + 50.69**2) / 12
    settings = (
        ("equilibrium", (1.531e9, 1.531e9, 0.1118e9)),
        (column_mass, (transverse, transverse, column_mass * radius**2 / 2)),
    )
    for number, ((mass, inertia), source_force) in enumerate(
            zip(settings, source_forces, strict=True), start=1):
        body = _body_data(HYDRO, number, mass, inertia, frequency,
                          1000, 9.81, .01)
        predicted = 2.25 * ramp[:, None] * (
            cosine[:, None] * body["excitation_real"]
            - sine[:, None] * body["excitation_imag"]
        )
        assert np.max(np.abs(source_force)) > 1e5
        error = np.max(np.abs(predicted - source_force))
        assert error < 1e-5, f"body {number} excitation gap: {error}"

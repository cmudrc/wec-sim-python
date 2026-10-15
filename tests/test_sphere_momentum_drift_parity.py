"""Full-duration public-WEC gate for the Sphere momentum-drift option."""

from hashlib import sha256
from pathlib import Path

import numpy as np

from wecsim import RegularCICWave, WEC


FIXTURE = Path(__file__).parent / "fixtures/matlab/sphere_momentum_drift"
HYDRO = Path(__file__).parent / "fixtures/matlab/sphere_mean_drift/sphere.h5"


def test_public_sphere_momentum_drift_against_local_matlab():
    source_file = FIXTURE / "source.npz"
    assert sha256(source_file.read_bytes()).hexdigest() == (
        "bc1a32a4b1e5f78b35396952754b80517825959553aff9898279adf9a76e768c"
    )
    with np.load(source_file, allow_pickle=False) as saved:
        body, wave, coefficients, mass = (
            saved[name].copy() for name in ("body", "wave", "coefficients", "mass")
        )
    assert (body.shape, wave.shape, coefficients.shape, mass.shape) == (
        (10_001, 43), (10_001, 2), (6, 3), (4,),
    )
    assert all(np.isfinite(values).all() for values in
               (body, wave, coefficients, mass))
    assert coefficients[0, 2] > 0  # The alternate mean-drift force is active.

    wec = WEC("Sphere momentum drift")
    sphere = wec.body("sphere", HYDRO, inertia=mass[1:],
                      mean_drift="momentum_conservation")
    for axis in ("surge", "heave", "pitch"):
        wec.coordinate(axis, sphere.move(axis))
    result = wec.run(RegularCICWave(0.1, 2), dt=0.01, end_time=100,
                     ramp_time=20, radiation_memory=10)
    response = result.raw
    assert result.case["bodies"][0]["mean_drift"] == "momentum_conservation"
    np.testing.assert_allclose(response.time, body[:, 0], rtol=0, atol=1e-8)
    np.testing.assert_allclose(response.wave_elevation, wave[:, 1],
                               rtol=0, atol=1e-12)

    # Reuse the published control-surface case's limits for this one-setting
    # source variant. No source motion or force is fed into the Python run.
    for axis, position_limit, speed_limit in (
        (0, 1e-4, 1e-4),       # surge
        (2, 7e-4, 2e-3),       # heave
        (4, 1e-6, 1e-6),       # pitch
    ):
        assert np.max(np.abs(response.body_position[:, 0, axis]
                             - body[:, 1 + axis])) < position_limit
        assert np.max(np.abs(response.body_velocity[:, 0, axis]
                             - body[:, 7 + axis])) < speed_limit
    np.testing.assert_allclose(response.body_position[:, 0, [1, 3, 5]],
                               body[:, [2, 4, 6]], rtol=0, atol=1e-10)
    np.testing.assert_allclose(response.body_velocity[:, 0, [1, 3, 5]],
                               body[:, [8, 10, 12]], rtol=0, atol=1e-10)

    amplitude = 0.05
    ramp = np.where(response.time < 20,
                    (1 - np.cos(np.pi * response.time / 20)) / 2, 1)
    phase = np.pi * response.time
    source_excitation = ramp[:, None] * (
        amplitude * (np.cos(phase[:, None]) * coefficients[:, 0]
                     - np.sin(phase[:, None]) * coefficients[:, 1])
        + amplitude**2 * coefficients[:, 2]
    )
    np.testing.assert_allclose(source_excitation, body[:, 19:25],
                               rtol=0, atol=2e-9)
    outputs = dict(response.extra_outputs)
    np.testing.assert_allclose(outputs["body1_excitation_force"],
                               body[:, 19:25], rtol=0, atol=2e-9)
    np.testing.assert_allclose(outputs["body1_mean_drift_force"],
                               amplitude**2 * ramp[:, None] * coefficients[:, 2],
                               rtol=0, atol=2e-9)

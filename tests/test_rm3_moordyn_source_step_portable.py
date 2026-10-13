"""Lock down the source-step sensitivity of published RM3 MoorDyn motion."""

from pathlib import Path

import numpy as np

from wecsim.environment.irregularWave import (
    jonswap_equal_energy_components, synthesize_irregular_response,
)


FIXTURE = (Path(__file__).parent / "fixtures/matlab"
           / "rm3_moordyn_source_step.npz")
HYDRO = Path(__file__).resolve().parents[1] / "examples/data/rm3.h5"


def _squared_gap_reduction(published, maxstep_only, refined):
    """Fraction of the published-to-refined squared gap removed by MaxStep."""
    original = refined - published
    residual = refined - maxstep_only
    assert np.sum(original**2) > 0
    return 1 - np.sum(residual**2) / np.sum(original**2)


def test_published_rm3_moordyn_surge_gap_is_source_maxstep_sensitive():
    """Compare three MATLAB runs on the same sea, without prescribing Python."""
    with np.load(FIXTURE, allow_pickle=False) as source:
        np.testing.assert_allclose(source["time"], np.arange(801) * .1,
                                   rtol=0, atol=1e-12)
        phase = source["phase"]
        wave = source["wave"]
        assert phase.shape == (3, 1000)
        assert wave.shape == (3, 801)
        np.testing.assert_array_equal(phase[0], phase[1])
        np.testing.assert_array_equal(phase[1], phase[2])
        np.testing.assert_allclose(wave[1], wave[0], rtol=0, atol=2e-14)
        np.testing.assert_allclose(wave[2], wave[0], rtol=0, atol=2e-14)

        for name, axis, original_minimum, reduction_minimum in (
            ("body1_position", 0, .1, .99),
            ("body2_position", 0, .1, .99),
            ("body1_position", 4, .001, .99),
            ("mooring_force", 4, 50_000, .99),
        ):
            history = source[name]
            assert history.shape == (3, 801, 6)
            published, maxstep_only, refined = history[:, :, axis]
            assert np.max(np.abs(refined - published)) > original_minimum
            reduction = _squared_gap_reduction(
                published, maxstep_only, refined,
            )
            assert reduction > reduction_minimum, (name, axis, reduction)

        # Shared joint pitch is the same on both bodies in every source run.
        np.testing.assert_allclose(
            source["body1_position"][:, :, 4],
            source["body2_position"][:, :, 4], rtol=0, atol=2e-15,
        )
        # The inference is channel-specific: spar heave does not show the
        # near-complete MaxStep explanation observed for surge and pitch.
        heave = source["body2_position"][:, :, 2]
        assert _squared_gap_reduction(*heave) < .6


def test_published_rm3_moordyn_incident_sea_from_pinned_phases():
    """Rebuild the 80 s sea independently of the saved body trajectories."""
    with np.load(FIXTURE, allow_pickle=False) as source:
        phase = source["phase"][0, :, None]
        elevation = source["wave"][0]
    components = jonswap_equal_energy_components(
        HYDRO, significant_height=2, peak_period=8,
        directions=np.array([0.]), spreading=np.array([1.]),
        phase=phase, discretization="traditional",
    )
    assert components.omega.shape == (1000,)
    wave = synthesize_irregular_response(
        HYDRO, components, dt=.1, end_time=80, ramp_time=0,
        rho=1025, g=9.81,
    )
    np.testing.assert_allclose(wave.time, np.arange(801) * .1,
                               rtol=0, atol=1e-12)
    np.testing.assert_allclose(wave.elevation, elevation,
                               rtol=0, atol=2e-14)

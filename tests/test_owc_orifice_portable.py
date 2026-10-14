"""Pinned OWC orifice law and pre-threshold physical motion in ordinary CI."""

from pathlib import Path

import numpy as np

from wecsim import OrificePTO, PMWave, WEC
from wecsim.environment.irregularWave import pm_equal_energy_components


FIXTURES = Path(__file__).parent / "fixtures/matlab"
HYDRO = FIXTURES / "owc_orifice.h5"
SOURCE = FIXTURES / "owc_orifice_prefix.npz"


def _limit(actual, expected, maximum, name):
    assert actual.shape == expected.shape, name
    assert np.isfinite(actual).all() and np.isfinite(expected).all(), name
    difference = float(np.max(np.abs(actual - expected)))
    assert difference < maximum, f"{name}: {difference:.8g} >= {maximum:.8g}"


def test_published_orifice_law_at_both_flow_signs_and_mach_regimes():
    with np.load(SOURCE, allow_pickle=False) as saved:
        parameters = saved["parameters"]
        samples = saved["orifice_samples"]
        np.testing.assert_array_equal(
            saved["orifice_sample_indices"], [0, 12, 263, 4330, 294352],
        )
    assert samples.shape == (5, 6)
    assert samples[1, 4] < 0 < samples[2, 4]
    np.testing.assert_array_equal(samples[:, 2], [0, 0, 0, 1, 1])
    model = OrificePTO(*parameters)
    piston_speed = samples[:, 4] / model.piston_area
    response = model.evaluate(piston_speed)
    _limit(response.force, samples[:, 1], 1e-8, "orifice reaction")
    _limit(response.pressure_drop / 1000, samples[:, 3], 1e-10,
           "orifice pressure in source kPa")
    _limit(response.flow_rate, samples[:, 4], 1e-13, "orifice flow")
    _limit(response.absorbed_power / 1000, samples[:, 5], 1e-10,
           "orifice power in source kW")
    np.testing.assert_array_equal(response.compressibility_flag, samples[:, 2])
    np.testing.assert_allclose(response.absorbed_power,
                               -response.force * piston_speed,
                               rtol=0, atol=1e-8)


def test_public_coupled_owc_motion_before_source_mach_flag():
    with np.load(SOURCE, allow_pickle=False) as saved:
        source = {name: saved[name] for name in
                  ("body", "flex", "wave", "components", "parameters")}
    body, flex, wave, components = (
        source[name] for name in ("body", "flex", "wave", "components")
    )
    assert body.shape == (1201, 25)
    assert flex.shape == (1201, 10)
    assert wave.shape == (1201, 2)
    assert components.shape == (500, 4)
    _limit(body[:, 0], np.arange(1201) * .005, 1e-10, "source time")
    _limit(wave[:, 0], body[:, 0], 1e-10, "source wave time")

    # The published input specifies substream 128. Generate its phases from
    # that setting; the saved MATLAB phases are comparison data only.
    phase = components[:, 3]
    incident = pm_equal_energy_components(
        HYDRO, significant_height=1, peak_period=4,
        directions=np.array([0.]), spreading=np.array([1.]),
        seed=128, phase_generator="matlab",
    )
    _limit(incident.phase[:, 0], phase, 1e-12, "published PM phases")
    _limit(incident.omega, components[:, 0], 1e-12, "PM frequency")
    _limit(incident.spectral_amplitude, components[:, 1], 1e-12,
           "PM spectral amplitude")
    _limit(incident.d_omega, components[:, 2], 1e-12, "PM bin width")
    wec = WEC("OWC physical prefix")
    owc = wec.body("OWC", HYDRO, inertia=(99.28, 11.04, 99.2))
    wec.floating_gbm(
        owc, orifice=OrificePTO(*source["parameters"]),
        heave_linear_damping=100, mode_linear_damping=100,
        heave_drag_cd=1.2, heave_drag_area=8,
        pitch_drag_cd=1.2, pitch_drag_area=8,
    )
    result = wec.run(
        PMWave(height=1, period=4, seed=128, phase_generator="matlab"),
        dt=.005, end_time=6, ramp_time=10, radiation_memory=15,
    )
    assert result.case["wave"]["seed"] == 128
    assert result.case["wave"]["phase_generator"] == "matlab"
    assert result.case["constraint"]["orifice_force_path"] == "coupled"
    _limit(result.wave_elevation, wave[:, 1], 1e-11, "incident sea")
    _limit(result.bodies["OWC"].position[:, 0], source["body"][:, 1],
           .012, "early surge")
    _limit(result.bodies["OWC"].position[:, 2], source["body"][:, 3],
           .025, "early heave")
    source_pitch = np.r_[0, np.cumsum(
        (body[1:, 11] + body[:-1, 11]) * .005 / 2,
    )]
    _limit(result.coordinates["pitch_unwrapped"].position, source_pitch,
           .010, "early pitch")
    _limit(result.flexible_modes["OWC"].position[:, 0], flex[:, 0],
           .0075, "early flexible mode")
    pto = result.ptos["orifice"]
    np.testing.assert_allclose(pto.absorbed_power, -pto.force * pto.velocity,
                               rtol=0, atol=1e-8)
    assert not np.any(dict(result.raw.extra_outputs)[
        "orifice_compressibility_flag"])

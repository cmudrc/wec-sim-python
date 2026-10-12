"""Short public hinge runs that preserve independently constructed seas."""

from pathlib import Path

import numpy as np
from scipy.io import savemat

from wecsim.dynamics.caseDynamics import run_case
from wecsim.dynamics.hingePitch import solve_hinged_pitch_from_excitation
from wecsim.environment.irregularWave import (
    imported_full_directional_components, pm_equal_energy_components,
    synthesize_full_directional_response, synthesize_multiple_irregular_response,
)


OSWEC = Path(__file__).resolve().parents[1] / "examples/data/oswec.h5"


def _hinge_case(wave):
    return {
        "simulation": {"dt": .05, "end_time": .2, "ramp_time": 0,
                       "radiation_memory": .1},
        "wave": wave,
        "bodies": [
            {"hydro_file": str(OSWEC), "mass": 127_000,
             "pitch_inertia": 1.85e6},
            {"name": "base", "nonhydro": True, "fixed": True,
             "center_gravity": [0, 0, -10.9], "mass": 999,
             "inertia": [1, 1, 1]},
        ],
        "constraint": {"kind": "fixed_hinge", "location": [0, 0, -10]},
        "pto": {"kind": "pitch", "location": [0, 0, -8.9],
                "damping": 120_000},
    }


def _assert_hinge_follows_incident(response, incident):
    solved = solve_hinged_pitch_from_excitation(
        OSWEC, incident.excitation_force, hinge_z=-8.9,
        body_mass=127_000, pitch_inertia=1.85e6,
        pto_damping=120_000, dt=.05, memory_time=.1,
    )
    np.testing.assert_allclose(response.wave_elevation, incident.elevation,
                               rtol=0, atol=1e-12)
    np.testing.assert_allclose(response.body_position[:, 0, 4], solved.angle,
                               rtol=0, atol=1e-11)
    np.testing.assert_allclose(response.pto_force, solved.pto_torque,
                               rtol=0, atol=1e-8)


def test_public_hinge_combines_two_pm_seas_before_solving_motion():
    seas = [
        {"height": 1., "period": 8., "direction": 0,
         "seed": 3, "frequency_count": 8},
        {"height": .5, "period": 8., "direction": 90,
         "seed": 4, "frequency_count": 8},
    ]
    response = run_case(_hinge_case({"type": "pm_multi", "seas": seas}))
    components = tuple(pm_equal_energy_components(
        OSWEC, significant_height=sea["height"],
        peak_period=sea["period"], directions=[sea["direction"]],
        spreading=[1], count=8, seed=sea["seed"],
    ) for sea in seas)
    incident = synthesize_multiple_irregular_response(
        OSWEC, components, dt=.05, end_time=.2, ramp_time=0,
    )
    _assert_hinge_follows_incident(response, incident)
    assert np.isfinite(response.body_position).all()
    assert len(response.time) == 5
    np.testing.assert_array_equal(response.body_velocity[:, 1], 0)


def test_public_hinge_imports_full_directional_spectrum(tmp_path):
    frequencies = np.array([.1, .15, .2])
    directions = np.array([0., 30., 60.])
    width = np.deg2rad(30)
    spectrum = tmp_path / "full_directional.mat"
    phase_file = tmp_path / "phase.csv"
    phases = np.zeros((3, 3))
    np.savetxt(phase_file, phases, delimiter=",")
    savemat(spectrum, {
        "frequencies": frequencies, "spectrum": np.array([.2, .4, .2]),
        "directions": directions, "spread": np.full((3, 3), 1 / (3 * width)),
    })
    response = run_case(_hinge_case({
        "type": "spectrumImportFullDir", "file": spectrum.name,
        "phase_file": phase_file.name,
    }), base_dir=tmp_path)
    components = imported_full_directional_components(
        OSWEC, spectrum, phase=phases,
    )
    incident = synthesize_full_directional_response(
        OSWEC, components, dt=.05, end_time=.2, ramp_time=0,
    )
    _assert_hinge_follows_incident(response, incident)
    assert set(response.auxiliary_files) == {spectrum, phase_file}
    assert np.isfinite(response.body_position).all()

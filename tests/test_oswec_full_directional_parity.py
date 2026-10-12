"""Separate physical heading quadrature from the pinned MATLAB force block."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim import FullDirectionalSpectrumWave, WEC, WorldPoint
from wecsim.dynamics.caseDynamics import run_case
from wecsim.environment.irregularWave import (
    imported_full_directional_components, synthesize_full_directional_response,
)


APPLICATIONS = os.environ.get("WEC_SIM_APPLICATIONS_DIR")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
pytestmark = pytest.mark.skipif(
    not (APPLICATIONS and REFERENCE
         and os.environ.get("WEC_SIM_REFERENCE_MODEL") == "OSWEC_FULL_DIR"),
    reason="paired MATLAB full-directional output not provided",
)


@pytest.fixture(scope="module")
def published():
    reference = Path(REFERENCE).resolve()
    applications = Path(APPLICATIONS).resolve()
    hydro = (applications / "_Common_Input_Files/OSWEC/hydroData/oswec.h5").resolve()
    spectrum = (applications / "Full_Directional_Waves/fullDirSpectrum.mat").resolve()
    phase_file = (reference / "OSWEC_FULL_DIR_phase.csv").resolve()
    phase = np.loadtxt(phase_file, delimiter=",")
    components = imported_full_directional_components(
        hydro, spectrum, seed=1, phase_generator="matlab",
    )
    assert np.max(np.abs(components.phase - phase)) < 1e-13
    bodies = tuple(np.loadtxt(
        reference / f"OSWEC_FULL_DIR_Full_Directional_Waves_body{number}.csv",
        delimiter=",",
    ) for number in (1, 2))
    pto = np.loadtxt(
        reference / "OSWEC_FULL_DIR_Full_Directional_Waves_pto1.csv",
        delimiter=",",
    )
    wave = np.loadtxt(reference / "OSWEC_FULL_DIR_wave.csv", delimiter=",")
    return hydro, spectrum, phase_file, components, bodies, pto, wave


def test_heading_width_explains_source_excitation_for_both_bodies(published):
    hydro, _, _, components, bodies, _, wave = published
    assert np.max(np.abs(components.d_theta - np.deg2rad(2))) < 1e-12
    for number, body in enumerate(bodies, start=1):
        source = synthesize_full_directional_response(
            hydro, components, dt=.05, end_time=400, ramp_time=100,
            body_number=number, excitation_interpolation="spline_frequency",
            force_quadrature="matlab_omitted",
        )
        integrated = synthesize_full_directional_response(
            hydro, components, dt=.05, end_time=400, ramp_time=100,
            body_number=number, excitation_interpolation="spline_frequency",
        )
        np.testing.assert_allclose(source.time, wave[:, 0], rtol=0, atol=1e-10)
        assert np.max(np.abs(source.elevation - wave[:, 1])) < 1e-11
        assert np.max(np.abs(integrated.elevation - wave[:, 1])) < 1e-11
        assert np.max(np.abs(source.excitation_force - body[:, 19:25])) < 1e-6
        # The pinned force block leaves out dTheta; integrated forcing must
        # carry the missing sqrt(dTheta) for this zero-mean-drift input.
        assert np.max(np.abs(
            integrated.excitation_force
            - body[:, 19:25] * np.sqrt(components.d_theta[0])
        )) < 1e-6


@pytest.mark.parametrize("scheme,limits", [
    ("implicit", (0.015, 0.004, 0.003,
                  0.025, 0.007, 0.005, 70, 50)),
    ("simulink_delay", (0.002, 0.0004, 0.0004,
                       0.003, 0.0007, 0.0006, 8, 5)),
])
def test_source_force_hinge_and_pto_against_matlab(published, scheme, limits):
    hydro, spectrum, _, _, bodies, pto, wave = published
    simulation = {"dt": .05, "end_time": 400, "ramp_time": 100,
                  "radiation_memory": 30}
    if scheme != "implicit":
        simulation["added_mass_scheme"] = scheme
    wave_config = FullDirectionalSpectrumWave(
        spectrum, seed=1, phase_generator="matlab",
        excitation_interpolation="spline_frequency",
        force_quadrature="matlab_omitted",
    )
    if scheme == "simulink_delay":
        configured = WEC("OSWEC full directional")
        flap = configured.body("flap", hydro, mass=127000,
                               inertia=(0, 1.85e6, 0), hydro_body=1)
        base = configured.body("base", hydro, mass=999,
                               inertia=(999, 999, 999), hydro_body=2)
        configured.fixed_hinge(
            flap, base, location=WorldPoint(0, 0, -10),
            pto_location=WorldPoint(0, 0, -8.9), damping=12000,
            added_mass_scheme=scheme,
        )
        response = configured.run(
            wave_config, dt=.05, end_time=400, ramp_time=100,
            radiation_memory=30,
        ).raw
    else:
        response = run_case({
            "simulation": simulation,
            "wave": wave_config.as_case(),
            "bodies": [
                {"hydro_file": str(hydro), "hydro_body": 1,
                 "mass": 127000, "pitch_inertia": 1.85e6},
                {"hydro_file": str(hydro), "hydro_body": 2, "fixed": True,
                 "mass": 999, "inertia": [999, 999, 999]},
            ],
            "constraint": {"kind": "fixed_hinge", "location": [0, 0, -10]},
            "pto": {"kind": "pitch", "location": [0, 0, -8.9],
                    "damping": 12000},
        })
    flap, base = bodies
    np.testing.assert_allclose(response.time, flap[:, 0], rtol=0, atol=1e-10)
    assert np.max(np.abs(response.wave_elevation - wave[:, 1])) < 1e-11
    for dof, pos_limit, vel_limit in (
        (0, limits[0], limits[3]),
        (2, limits[1], limits[4]),
        (4, limits[2], limits[5]),
    ):
        assert np.max(np.abs(response.body_position[:, 0, dof]
                             - flap[:, 1 + dof])) < pos_limit
        assert np.max(np.abs(response.body_velocity[:, 0, dof]
                             - flap[:, 7 + dof])) < vel_limit
    for dof in (1, 3, 5):
        np.testing.assert_allclose(response.body_position[:, 0, dof],
                                   flap[:, 1 + dof], rtol=0, atol=1e-10)
        np.testing.assert_allclose(response.body_velocity[:, 0, dof],
                                   flap[:, 7 + dof], rtol=0, atol=1e-10)
    np.testing.assert_allclose(response.body_position[:, 1], base[:, 1:7],
                               rtol=0, atol=1e-10)
    np.testing.assert_allclose(response.body_velocity[:, 1], base[:, 7:13],
                               rtol=0, atol=1e-10)
    assert np.max(np.abs(response.pto_force - pto[:, 17])) < limits[6]
    power = response.pto_force * response.body_velocity[:, 0, 4]
    assert np.max(np.abs(power - pto[:, 23])) < limits[7]
    assert len(response.hydro_files) == 2
    assert len(response.auxiliary_files) == 1

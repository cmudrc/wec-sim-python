"""Pair the published two-sea OSWEC forces, motion, and PTO with MATLAB."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim import MultiPMWave, PMWave, WEC, WorldPoint
from wecsim.caseDynamics import run_case
from wecsim.irregularWave import (
    IrregularComponents, synthesize_multiple_irregular_response,
)


APPLICATIONS = os.environ.get("WEC_SIM_APPLICATIONS_DIR")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
pytestmark = pytest.mark.skipif(
    not (APPLICATIONS and REFERENCE
         and os.environ.get("WEC_SIM_REFERENCE_MODEL") == "OSWEC_MULTI_WAVE"),
    reason="paired MATLAB multiple-wave output not provided",
)


@pytest.fixture(scope="module")
def published():
    reference = Path(REFERENCE)
    hydro = (Path(APPLICATIONS) / "_Common_Input_Files/OSWEC/hydroData/oswec.h5").resolve()
    components = []
    for number, direction in ((1, 0), (2, 90)):
        values = np.loadtxt(
            reference / f"OSWEC_MULTI_WAVE_wave{number}_components.csv",
            delimiter=",",
        )
        components.append(IrregularComponents(
            omega=values[:, 0], spectral_amplitude=values[:, 1],
            d_omega=values[:, 2], directions=np.array([direction]),
            spreading=np.array([1.0]), phase=values[:, 3, None],
        ))
    bodies = tuple(np.loadtxt(
        reference / f"OSWEC_MULTI_WAVE_Multiple_Wave_Spectra_body{number}.csv",
        delimiter=",",
    ) for number in (1, 2))
    pto = np.loadtxt(
        reference / "OSWEC_MULTI_WAVE_Multiple_Wave_Spectra_pto1.csv",
        delimiter=",",
    )
    return hydro, components, bodies, pto


def test_combined_elevation_and_both_body_excitation_forces(published):
    hydro, components, bodies, _ = published
    expected_elevation = sum(np.loadtxt(
        Path(REFERENCE) / f"OSWEC_MULTI_WAVE_wave{number}_elevation.csv",
        delimiter=",",
    )[:, 1] for number in (1, 2))
    for number, source in enumerate(bodies, start=1):
        incident = synthesize_multiple_irregular_response(
            hydro, components, dt=0.1, end_time=100, ramp_time=0.1,
            body_number=number, excitation_interpolation="spline",
        )
        np.testing.assert_allclose(incident.time, source[:, 0],
                                   rtol=0, atol=1e-10)
        assert np.max(np.abs(incident.elevation - expected_elevation)) < 1e-11
        assert np.max(np.abs(incident.excitation_force - source[:, 19:25])) < 1e-6


@pytest.mark.parametrize("scheme,limits", [
    ("implicit", (0.11, 0.012, 0.023, 0.32, 0.033, 0.065, 8000, 2600)),
    ("simulink_delay", (0.012, 0.0014, 0.0025,
                       0.04, 0.004, 0.008, 900, 350)),
])
def test_multiple_wave_hinge_and_pto_against_matlab(published, tmp_path,
                                                    scheme, limits):
    hydro, components, bodies, pto = published
    seas = []
    for number, height, direction in ((1, 2, 0), (2, 1, 90)):
        phase_file = tmp_path / f"phase{number}.csv"
        np.savetxt(phase_file, components[number - 1].phase, delimiter=",")
        seas.append({"height": height, "period": 3, "direction": direction,
                     "phase_file": phase_file.name})
    simulation = {"dt": 0.1, "end_time": 100, "ramp_time": 0.1,
                  "radiation_memory": 40}
    if scheme != "implicit":
        simulation["added_mass_scheme"] = scheme
    response = run_case({
        "simulation": simulation,
        "wave": {"type": "pm_multi", "seas": seas,
                 "excitation_interpolation": "spline"},
        "bodies": [
            {"hydro_file": str(hydro), "hydro_body": 1,
             "mass": 12700, "pitch_inertia": 1.85e6},
            {"hydro_file": str(hydro), "hydro_body": 2,
             "fixed": True, "mass": 999, "inertia": [999, 999, 999]},
        ],
        "constraint": {"kind": "fixed_hinge", "location": [0, 0, -10]},
        "pto": {"kind": "pitch", "location": [0, 0, -8.9],
                "damping": 120000},
    }, base_dir=tmp_path)
    flap, base = bodies
    np.testing.assert_allclose(response.time, flap[:, 0], rtol=0, atol=1e-10)
    np.testing.assert_allclose(response.time, pto[:, 0], rtol=0, atol=1e-10)
    for dof, position_limit, velocity_limit in (
        (0, limits[0], limits[3]),
        (2, limits[1], limits[4]),
        (4, limits[2], limits[5]),
    ):
        assert np.max(np.abs(
            response.body_position[:, 0, dof] - flap[:, 1 + dof],
        )) < position_limit
        assert np.max(np.abs(
            response.body_velocity[:, 0, dof] - flap[:, 7 + dof],
        )) < velocity_limit
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
    source_power = pto[:, 23]
    calculated_power = response.pto_force * response.body_velocity[:, 0, 4]
    assert np.max(np.abs(calculated_power - source_power)) < limits[7]
    assert len(response.hydro_files) == 2
    assert len(response.auxiliary_files) == 2
    if scheme == "simulink_delay":
        wec = WEC("Published two-sea OSWEC")
        flap_body = wec.body("flap", hydro, mass=12700,
                             inertia=(1.85e6,) * 3)
        base_body = wec.body("base", hydro, mass=999,
                             inertia=(999,) * 3)
        wec.fixed_hinge(
            flap_body, base_body, location=WorldPoint(0, 0, -10),
            pto_location=WorldPoint(0, 0, -8.9), damping=120000,
            added_mass_scheme="simulink_delay",
        )
        public_wave = MultiPMWave(tuple(
            PMWave(sea["height"], sea["period"],
                   direction=sea["direction"], phase_file=sea["phase_file"])
            for sea in seas
        ), excitation_interpolation="spline")
        public = wec.run(public_wave, dt=.1, end_time=100,
                         ramp_time=.1, radiation_memory=40,
                         base_dir=tmp_path)
        np.testing.assert_allclose(public.time, response.time,
                                   rtol=0, atol=1e-10)
        np.testing.assert_allclose(public.raw.body_position,
                                   response.body_position, rtol=0, atol=1e-10)
        np.testing.assert_allclose(public.raw.body_velocity,
                                   response.body_velocity, rtol=0, atol=1e-10)
        np.testing.assert_allclose(public.ptos["hinge"].force,
                                   response.pto_force, rtol=0, atol=1e-8)
        np.testing.assert_allclose(public.wave_elevation,
                                   response.wave_elevation, rtol=0, atol=1e-10)

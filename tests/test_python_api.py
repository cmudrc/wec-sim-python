"""Exercise the public Python WEC builder without JSON or CSV files."""

import json
from pathlib import Path

import numpy as np
import pytest
from scipy.io import loadmat

from examples.configurable_rm3_pto import HYDRO, build_wec
from wecsim.caseDynamics import run_case
from wecsim import (FullDirectionalSpectrumWave, ImportedElevationWave, ImportedSpectrumWave, JONSWAPWave, LatchingControl, NoWave, PMWave,
                    RegularWave, WEC, WorldPoint)
from wecsim.irregularWave import (imported_spectrum_components,
                                  jonswap_equal_energy_components,
                                  pm_equal_energy_components,
                                  synthesize_irregular_response)


ROOT = Path(__file__).resolve().parents[1]
JSON_EXAMPLE = ROOT / "examples/configurable_rm3_pto.json"


def test_full_directional_wave_keeps_physical_defaults():
    wave = FullDirectionalSpectrumWave("fullDirSpectrum.mat", seed=7)
    assert wave.as_case() == {
        "type": "spectrumImportFullDir", "file": "fullDirSpectrum.mat",
        "seed": 7, "excitation_interpolation": "linear",
        "force_quadrature": "integrated",
    }
    with pytest.raises(ValueError, match="either seed or phase_file"):
        FullDirectionalSpectrumWave("fullDirSpectrum.mat", seed=7,
                                    phase_file="phase.csv").as_case()
    with pytest.raises(ValueError, match="substream seed"):
        FullDirectionalSpectrumWave("fullDirSpectrum.mat",
                                    phase_generator="matlab").as_case()


def test_python_builder_matches_equivalent_case():
    wec = build_wec()
    wave = RegularWave(height=2.5, period=8)
    result = wec.run(wave, dt=0.1, end_time=0.2, ramp_time=5)
    case = json.loads(JSON_EXAMPLE.read_text(encoding="utf-8"))
    case["simulation"]["end_time"] = 0.2
    reference = run_case(case, base_dir=JSON_EXAMPLE.parent)
    np.testing.assert_allclose(result.time, reference.time)
    np.testing.assert_allclose(result.raw.body_position, reference.body_position)
    np.testing.assert_allclose(
        result.ptos["main"].force,
        dict(reference.extra_outputs)["pto_main_force"],
    )
    np.testing.assert_allclose(
        result.coordinates["shared_pitch"].position,
        dict(reference.extra_outputs)["coordinate_shared_pitch_position"],
    )
    assert result.bodies["float"].position.shape == (3, 6)
    assert result.ptos["main"].absorbed_power.shape == (3,)
    case_mapping = wec.to_case(wave, dt=0.1, end_time=0.2, ramp_time=5)
    assert "coordinates" in case_mapping["constraint"]
    assert "coordinate_map" not in case_mapping["bodies"][0]
    json.dumps(case_mapping)


def test_public_fixed_hinge_reports_named_motion_and_dissipated_power():
    hydro = (ROOT / "tests/test_objects/test_bodyclass/testData/hydroData/oswec.h5")
    wec = WEC("Hinged OSWEC")
    flap = wec.body("flap", hydro, mass=127_000,
                    inertia=(0, 1.85e6, 0))
    base = wec.fixed_body("base", center_gravity=(0, 0, -10.9),
                          mass=999, inertia=(1, 1, 1))
    wec.fixed_hinge(flap, base, damping=12_000, stiffness=500,
                    equilibrium_angle=0.1)
    result = wec.run(RegularWave(2.5, 8), dt=0.1, end_time=1,
                     ramp_time=0)
    assert result.case["constraint"]["kind"] == "fixed_hinge"
    assert result.case["pto"]["equilibrium_position"] == 0.1
    np.testing.assert_array_equal(result.coordinates["pitch"].position,
                                  result.ptos["hinge"].stroke)
    np.testing.assert_array_equal(result.coordinates["pitch"].velocity,
                                  result.ptos["hinge"].velocity)
    np.testing.assert_allclose(
        result.ptos["hinge"].force,
        -12_000 * result.ptos["hinge"].velocity
        - 500 * (result.ptos["hinge"].stroke - 0.1),
        rtol=0, atol=1e-10,
    )
    np.testing.assert_allclose(
        result.ptos["hinge"].absorbed_power,
        12_000 * result.ptos["hinge"].velocity**2,
        rtol=0, atol=1e-10,
    )
    np.testing.assert_array_equal(result.bodies["base"].position[:, 2],
                                  np.full(11, -10.9))


def test_python_builder_supports_ground_anchor_and_named_initial_state():
    wec = WEC("Anchored float")
    body = wec.body("float", HYDRO)
    wec.coordinate("heave", body.move("heave"))
    wec.pto("spring", body.at(0, 0, 0), WorldPoint(0, 0, 10),
            stiffness=1000)
    result = wec.run(
        NoWave(), dt=0.1, end_time=0.2, radiation_memory=0.2,
        initial_coordinate={"heave": 1},
    )
    np.testing.assert_allclose(result.ptos["spring"].stroke[0], -1)
    np.testing.assert_allclose(result.ptos["spring"].force[0], 1000)
    assert np.isfinite(result.bodies["float"].position).all()


def test_python_builder_runs_jonswap_sea():
    wec = WEC("JONSWAP float")
    body = wec.body("float", HYDRO)
    wec.coordinate("heave", body.move("heave"))
    wave = JONSWAPWave(2.5, 8, seed=1, gamma=3.3)
    assert wec.to_case(wave, dt=.1, end_time=.2)["wave"]["gamma"] == 3.3
    result = wec.run(wave, dt=.1, end_time=.2, ramp_time=1,
                     radiation_memory=.2)
    assert result.wave_elevation.shape == (3,)
    assert np.isfinite(result.bodies["float"].position).all()


@pytest.mark.parametrize("wave_type,builder", [
    (PMWave, pm_equal_energy_components),
    (JONSWAPWave, jonswap_equal_energy_components),
])
def test_python_builder_uses_custom_irregular_frequency_range(wave_type, builder):
    wec = WEC("Narrow-band float")
    body = wec.body("float", HYDRO)
    wec.coordinate("heave", body.move("heave"))
    wave = wave_type(2.5, 8, seed=7, frequency_count=32,
                     frequency_range=(0.5, 1.5))
    case = wec.to_case(wave, dt=0.1, end_time=0.2, ramp_time=1)
    assert case["wave"]["frequency_range"] == [0.5, 1.5]
    assert "water_depth" not in case["wave"]
    response = wec.run(wave, dt=0.1, end_time=0.2,
                       ramp_time=1, radiation_memory=0.2)
    components = builder(
        HYDRO, significant_height=2.5, peak_period=8,
        directions=[0], spreading=[1], count=32, seed=7,
        frequency_range=(0.5, 1.5),
    )
    assert 0.5 <= components.omega.min() < components.omega.max() <= 1.5
    expected = synthesize_irregular_response(
        HYDRO, components, dt=0.1, end_time=0.2, ramp_time=1,
    )
    np.testing.assert_allclose(response.wave_elevation, expected.elevation,
                               rtol=0, atol=1e-12)
    assert np.isfinite(response.bodies["float"].position).all()


def test_floating_joint_accepts_traditional_jonswap_from_public_wave():
    wec = WEC("RM3 traditional sea")
    float_body = wec.body("float", HYDRO, inertia=(0, 21_306_090.66, 0))
    spar = wec.body("spar", HYDRO, inertia=(0, 94_407_091.24, 0))
    wec.floating_joint(float_body, spar, damping=1_200_000,
                       radiation_method="convolution")
    wave = JONSWAPWave(2, 8, seed=1, phase_generator="matlab",
                       discretization="traditional", frequency_count=64)
    result = wec.run(wave, dt=0.1, end_time=0.3,
                     ramp_time=0, radiation_memory=0.2)
    assert result.case["wave"]["discretization"] == "traditional"
    assert result.case["wave"]["frequency_count"] == 64
    components = jonswap_equal_energy_components(
        HYDRO, significant_height=2, peak_period=8,
        directions=[0], spreading=[1], count=64, seed=1,
        phase_generator="matlab", discretization="traditional",
    )
    expected = synthesize_irregular_response(
        HYDRO, components, dt=0.1, end_time=0.3, ramp_time=0,
    )
    np.testing.assert_allclose(result.wave_elevation, expected.elevation,
                               rtol=0, atol=1e-12)
    assert np.isfinite(result.bodies["float"].position).all()
    assert np.max(np.abs(result.bodies["float"].position[:, 2])) > 0
    assert JONSWAPWave(2, 8, discretization="traditional").as_case()[
        "frequency_count"] == 1000


def test_hydrodynamic_frequency_range_clamps_to_bem_limits():
    settings = dict(significant_height=2.5, peak_period=8,
                    directions=[0], spreading=[1], count=32, seed=7)
    default = pm_equal_energy_components(HYDRO, **settings)
    clamped = pm_equal_energy_components(
        HYDRO, frequency_range=(0.001, 10), **settings,
    )
    np.testing.assert_array_equal(clamped.omega, default.omega)
    np.testing.assert_array_equal(clamped.d_omega, default.d_omega)


def test_hydrodynamic_wave_rejects_unimplemented_depth_override():
    wec = WEC("Depth override")
    body = wec.body("float", HYDRO)
    wec.coordinate("heave", body.move("heave"))
    with pytest.raises(ValueError, match="fixed Morison body"):
        wec.run(PMWave(2.5, 8, seed=7, water_depth=25),
                dt=0.1, end_time=0.2, ramp_time=1,
                radiation_memory=0.2)


def test_python_builder_runs_imported_spectrum_with_relative_mat_path():
    source = ("tests/test_objects/test_waveclass/testData/"
              "spectrumImport_1_test/spectrumData.mat")
    wec = build_wec()
    wave = ImportedSpectrumWave(source)
    response = wec.run(
        wave, dt=0.1, end_time=0.2, ramp_time=1,
        radiation_memory=0.2, base_dir=ROOT,
    )
    components = imported_spectrum_components(HYDRO, ROOT / source)
    expected = synthesize_irregular_response(
        HYDRO, components, dt=0.1, end_time=0.2, ramp_time=1,
    )
    np.testing.assert_allclose(response.wave_elevation, expected.elevation,
                               rtol=0, atol=1e-12)
    assert response.raw.auxiliary_files == ((ROOT / source).resolve(),)
    assert np.isfinite(response.bodies["float"].position).all()
    assert np.isfinite(response.bodies["spar"].position).all()
    assert np.isfinite(response.ptos["main"].force).all()


def test_python_builder_runs_imported_elevation_with_relative_mat_path():
    source = ("tests/test_objects/test_waveclass/testData/"
              "etaImport_1_test/etaData.mat")
    wec = build_wec()
    response = wec.run(
        ImportedElevationWave(source), dt=0.1, end_time=0.2,
        ramp_time=1, radiation_memory=0.2, base_dir=ROOT,
    )
    raw = loadmat(ROOT / source)["etaData"]
    expected = np.interp(response.time, raw[:, 0], raw[:, 1])
    ramp = (1 - np.cos(np.pi * response.time)) / 2
    np.testing.assert_allclose(response.wave_elevation, expected * ramp,
                               rtol=0, atol=1e-12)
    assert response.raw.auxiliary_files == ((ROOT / source).resolve(),)
    assert np.isfinite(response.bodies["float"].position).all()
    assert np.isfinite(response.bodies["spar"].position).all()
    assert np.isfinite(response.ptos["main"].force).all()
    assert dict(response.raw.extra_outputs)["body1_excitation_force"].shape == (3, 6)


def test_python_builder_configures_paired_floating_joint_and_pto():
    wec = WEC("RM3 floating joint")
    float_body = wec.body("float", HYDRO, inertia=(0, 21_306_090.66, 0))
    spar = wec.body("spar", HYDRO, inertia=(0, 94_407_091.24, 0))
    wec.floating_joint(
        float_body, spar, damping=1_200_000, stiffness=100_000,
        equilibrium_position=0.1, pto_name="main",
    )
    result = wec.run(RegularWave(2.5, 8), dt=0.1, end_time=0.2,
                     ramp_time=1)
    assert set(result.coordinates) == {"surge", "float_heave",
                                       "spar_heave", "pitch"}
    pto = result.ptos["main"]
    np.testing.assert_allclose(
        pto.force,
        -1_200_000 * pto.velocity - 100_000 * (pto.stroke - 0.1),
    )
    np.testing.assert_allclose(pto.absorbed_power,
                               1_200_000 * pto.velocity**2)
    np.testing.assert_allclose(
        pto.stroke,
        result.coordinates["float_heave"].position
        - result.coordinates["spar_heave"].position,
    )


def test_python_builder_rejects_foreign_body_attachment():
    first = WEC("first")
    second = WEC("second")
    first.body("body", HYDRO)
    other = second.body("body", HYDRO)
    with pytest.raises(ValueError, match="must belong to this WEC"):
        first.pto("bad", other.at(0, 0, 0), WorldPoint(0, 0, 2),
                  damping=1)


def test_python_builder_exposes_latching_pto_settings_and_force():
    wec = WEC("Latching float")
    body = wec.body("float", HYDRO)
    wec.coordinate("heave", body.move("heave"))
    control = LatchingControl(gain=100, latch_damping=10_000,
                             latch_time=0.2)
    wec.pto("latch", body.at(0, 0, 0), WorldPoint(0, 0, 10),
            axis=(0, 0, 1), control=control)
    wave = RegularWave(2.5, 8)
    case = wec.to_case(wave, dt=0.01, end_time=0.1,
                       initial_speed={"heave": 0.1})
    assert case["ptos"][0]["control"] == {
        "kind": "latching", "latch_time": 0.2, "latch_damping": 10_000,
        "minimum_normal_time": 0.2,
    }
    result = wec.run(wave, dt=0.01, end_time=0.1, ramp_time=1,
                     initial_speed={"heave": 0.1})
    np.testing.assert_allclose(result.ptos["latch"].force[0],
                               -100 * result.ptos["latch"].velocity[0])
    assert np.isfinite(result.bodies["float"].position).all()

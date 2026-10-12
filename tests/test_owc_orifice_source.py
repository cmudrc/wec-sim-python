"""Pair the published OWC orifice, force paths, and early motion with MATLAB."""

import os
from pathlib import Path
from xml.etree import ElementTree
from zipfile import ZipFile

import numpy as np
import pytest
from scipy.io import loadmat
from scipy.signal import fftconvolve

from wecsim import OrificePTO, PMWave, WEC
from wecsim.hydrodynamics.bodyClass import BodyClass
from wecsim.environment.irregularWave import (
    pm_equal_energy_components, synthesize_irregular_response,
)


REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
APPLICATIONS = os.environ.get("WEC_SIM_APPLICATIONS_DIR")
APPLICATION_SOURCE = os.environ.get("WEC_SIM_APPLICATIONS_SOURCE_DIR", APPLICATIONS)
pytestmark = pytest.mark.skipif(
    not REFERENCE, reason="fresh MATLAB OWC orifice output absent",
)


def _max_error(actual, expected, limit, label):
    assert actual.shape == expected.shape, label
    assert np.isfinite(actual).all() and np.isfinite(expected).all(), label
    error = float(np.max(np.abs(actual - expected)))
    assert error < limit, f"{label}: {error:.6g} exceeds {limit}"


def _load(name):
    return np.loadtxt(Path(REFERENCE) / f"OWC_ORIFICE_{name}.csv",
                      delimiter=",", ndmin=2)


def test_published_orifice_force_on_matlab_piston_motion():
    reference = Path(REFERENCE)
    logged = _load("orifice")
    parameters = _load("parameters").ravel()
    flex = loadmat(reference / "OWC_ORIFICE_Flex_out.mat",
                   simplify_cells=True)["Flex_out"]
    time = np.asarray(flex["time"])
    values = np.asarray(flex["signals"]["values"])
    assert logged.ndim == 2 and logged.shape[1] == 6
    assert values.ndim == 2 and values.shape[1] == 10
    assert len(time) == len(values) and np.all(np.diff(time) > 0)
    assert np.all(np.diff(logged[:, 0]) > 0)
    model = OrificePTO(*parameters)
    # The orifice block logs every ode23t internal evaluation (418,828 rows),
    # while Flex_out is sampled at 0.005 s (26,001 rows). The logged flow gives
    # the exact internal piston speed. Its interpolated value agrees with the
    # independent flexible-mode speed on the output grid.
    piston_speed = logged[:, 4] / model.piston_area
    sampled_speed = np.interp(time, logged[:, 0], piston_speed)
    _max_error(sampled_speed, values[:, 1], 1e-10,
               "orifice input versus flexible-mode speed")
    result = model.evaluate(piston_speed)

    # MATLAB logs [force N, Mach flag, pressure kPa, flow m^3/s, power kW].
    _max_error(result.force, logged[:, 1], 1e-8, "orifice reaction")
    np.testing.assert_array_equal(result.compressibility_flag, logged[:, 2])
    _max_error(result.pressure_drop / 1000, logged[:, 3], 1e-10,
               "orifice pressure")
    _max_error(result.flow_rate, logged[:, 4], 1e-13, "orifice flow")
    _max_error(result.absorbed_power / 1000, logged[:, 5], 1e-10,
               "orifice power")
    _max_error(result.absorbed_power, -result.force * piston_speed,
               1e-8, "passive orifice energy")


@pytest.mark.skipif(not APPLICATION_SOURCE,
                    reason="published OWC Simulink model absent")
def test_owc_source_piston_omitted_from_flexible_state_integration():
    model = (Path(APPLICATION_SOURCE) / "OWC/OrificeModel/OWC_GBM.slx")
    with ZipFile(model) as archive:
        system = ElementTree.fromstring(archive.read(
            "simulink/systems/system_831_5393.xml",
        ))
    edges = set()
    for line in system.findall(".//Line"):
        ports = [(item.get("Name"), item.text)
                 for item in line.iter("P")]
        source = next(value for name, value in ports if name == "Src")
        edges.update((source, value) for name, value in ports
                     if name == "Dst")
    # The state-space input is mass^-1 times hydrodynamic force. The piston
    # force joins the separate reported-acceleration path instead.
    assert ("831:5394#out:1", "831:5408#in:2") in edges
    assert ("831:5408#out:1", "831:5404#in:2") in edges
    assert ("831:5404#out:1", "831:5401#in:1") in edges
    assert ("831:5402#out:1", "831:5395#in:4") in edges
    assert ("831:5395#out:1", "831:5411#in:2") in edges

    flexible = loadmat(
        Path(REFERENCE) / "OWC_ORIFICE_Flex_out.mat",
        simplify_cells=True,
    )["Flex_out"]
    time = np.asarray(flexible["time"]).ravel()
    values = np.asarray(flexible["signals"]["values"])
    mass = _load("flexible_effective_mass").item()
    velocity_derivative = np.gradient(values[:, 1], time)
    early = (time > 0) & (time < 6)
    hydro_only_error = velocity_derivative[early] - values[early, 3] / mass
    reported_error = velocity_derivative[early] - values[early, 2]
    assert np.sqrt(np.mean(hydro_only_error**2)) < 0.0021
    assert np.max(np.abs(hydro_only_error)) < 0.021
    assert np.sqrt(np.mean(reported_error**2)) > 0.25
    assert np.max(np.abs(reported_error)) > 0.7


@pytest.mark.skipif(not APPLICATIONS,
                    reason="BEMIO-generated OWC HDF5 absent")
def test_published_owc_wave_and_seventh_excitation_channel():
    hydro = (Path(APPLICATIONS)
             / "OWC/OrificeModel/hydroData/test17a_clean.h5")
    components = _load("components")
    wave = _load("wave")
    rigid = _load("OrificeModel_body1")
    flex = loadmat(Path(REFERENCE) / "OWC_ORIFICE_Flex_out.mat",
                   simplify_cells=True)["Flex_out"]
    flexible = np.asarray(flex["signals"]["values"])
    assert components.shape == (500, 4)
    assert wave.shape == (26001, 2)
    assert rigid.shape == (26001, 25)
    assert flexible.shape == (26001, 10)

    realized = pm_equal_energy_components(
        hydro, significant_height=1, peak_period=4,
        directions=np.array([0.0]), spreading=np.array([1.0]),
        phase=components[:, 3, None],
    )
    _max_error(realized.omega, components[:, 0], 1e-12, "PM frequency")
    _max_error(realized.spectral_amplitude, components[:, 1], 1e-12,
               "PM spectral amplitude")
    _max_error(realized.d_omega, components[:, 2], 1e-12, "PM bin width")
    result = synthesize_irregular_response(
        hydro, realized, dt=0.005, end_time=130, ramp_time=10,
    )
    assert result.excitation_force.shape == (26001, 7)
    _max_error(result.time, wave[:, 0], 1e-10, "wave time")
    _max_error(result.elevation, wave[:, 1], 1e-11, "wave elevation")
    _max_error(result.excitation_force[:, :6], rigid[:, 19:25], 1e-7,
               "rigid excitation")
    _max_error(result.excitation_force[:, 6], flexible[:, 4], 1e-7,
               "flexible excitation")

    body = BodyClass(str(hydro))
    body.bodyNumber = body.bodyTotal = 1
    body.readH5file()
    body.mass = "equilibrium"
    body.linearDamping = np.diag([0, 0, 100, 0, 0, 0, 100])
    memory = np.arange(3001) * 0.005
    body.hydroForcePre(
        components[:, 0], [0], len(memory), memory, len(components),
        0.005, 1000, 9.81, "irregular", np.zeros((2, 1)), 1, 1, 0, 0, 0,
    )
    force = body.hydroForce
    for key, filename in (
        ("mass_ff", "flexible_effective_mass"),
        ("stiffness", "flexible_stiffness"),
        ("damping", "flexible_damping"),
    ):
        _max_error(np.asarray(force["gbm"][key]), _load(filename), 1e-9,
                   f"flexible {key}")
    _max_error(np.asarray(force["fAddedMass"]), _load("added_mass"),
               1e-9, "added mass")
    _max_error(np.asarray(force["linearHydroRestCoef"]),
               _load("hydrostatic_stiffness"), 1e-9,
               "hydrostatic stiffness")

    # The source flexible state equation applies the orifice reaction in
    # addition to its logged hydrodynamic force balance.
    total = flexible[:, 3]
    hydro_terms = flexible[:, 4:]
    _max_error(total, hydro_terms[:, 0] - hydro_terms[:, 1:].sum(axis=1),
               1e-8, "source flexible hydrodynamic force balance")
    source_orifice = _load("orifice")
    reaction = np.interp(wave[:, 0], source_orifice[:, 0], source_orifice[:, 1])
    mass = float(np.asarray(force["gbm"]["mass_ff"])[0, 0])
    _max_error(mass * flexible[:, 2], total + reaction, 1e-6,
               "source flexible acceleration and orifice balance")

    # Reconstruct radiation from all seven source speeds.
    # The finite-memory convolution is trapezoidal at its first and last
    # samples, matching the source's continuous convolution block.
    velocities = np.column_stack((rigid[:, 7:13], flexible[:, 1]))
    kernel = np.asarray(force["irkb"])
    radiation = 0.005 * np.column_stack(
        [sum(
            fftconvolve(kernel[:, axis, channel], velocities[:, channel])
            [:len(wave)] for channel in range(7)
        ) for axis in range(7)]
    )
    radiation -= 0.0025 * velocities @ kernel[0].T
    last_lag = len(kernel) - 1
    radiation[last_lag:] -= 0.0025 * velocities[:-last_lag] @ kernel[-1].T
    _max_error(radiation[:, 6], flexible[:, 5], 0.01,
               "source seventh radiation force")

    rigid_forces = _load("body1_forces")
    assert rigid_forces.shape == (26001, 37)
    rigid_radiation, rigid_added, rigid_restoring, rigid_viscous, \
        rigid_linear, rigid_acceleration = (
            rigid_forces[:, 1 + 6*i:1 + 6*(i+1)] for i in range(6)
        )
    _max_error(radiation[:, :6], rigid_radiation, 0.25,
               "source rigid radiation force")
    displacements = np.column_stack((rigid[:, 1:7], flexible[:, 0]))
    _max_error(
        displacements @ np.asarray(force["linearHydroRestCoef"])[:6].T,
        rigid_restoring, 1e-7, "source rigid hydrostatic force",
    )
    _max_error(
        velocities @ np.asarray(force["linearDamping"])[:6].T,
        rigid_linear, 1e-8, "source rigid linear damping",
    )
    _max_error(
        rigid[:, 19:25] - rigid_radiation - rigid_added
        - rigid_restoring - rigid_viscous - rigid_linear,
        rigid[:, 13:19], 1e-7, "source rigid hydrodynamic force sum",
    )
    assert np.isfinite(rigid_acceleration).all()
    rigid_mass = _load("body1_mass").ravel()
    shifted_mass = (rigid_mass[0]
                    + 2 * np.trace(np.asarray(force["fAddedMass"])[:3, :3]))
    rigid_reaction = np.interp(
        rigid[:, 0], source_orifice[:, 0], source_orifice[:, 1],
    )
    _max_error(shifted_mass * rigid_acceleration[:, 0],
               rigid[:, 13], 1e-7, "source rigid surge balance")
    _max_error(shifted_mass * rigid_acceleration[:, 2] + rigid_reaction,
               rigid[:, 15], 1e-7, "source rigid heave and piston balance")
    _max_error(rigid_mass[2] * rigid_acceleration[:, 4],
               rigid[:, 17], 1e-8, "source rigid pitch balance")


@pytest.mark.skipif(not APPLICATIONS,
                    reason="BEMIO-generated OWC HDF5 absent")
def test_owc_coupled_motion_before_source_air_flag(tmp_path):
    hydro = (Path(APPLICATIONS)
             / "OWC/OrificeModel/hydroData/test17a_clean.h5")
    source = _load("OrificeModel_body1")
    flexible = loadmat(
        Path(REFERENCE) / "OWC_ORIFICE_Flex_out.mat",
        simplify_cells=True,
    )["Flex_out"]["signals"]["values"]
    components = _load("components")
    phase_file = tmp_path / "owc_phase.csv"
    np.savetxt(phase_file, components[:, 3], delimiter=",")
    model = WEC("OWC")
    body = model.body("OWC", hydro, inertia=(99.28, 11.04, 99.2))
    orifice = OrificePTO(*_load("parameters").ravel())
    model.floating_gbm(
        body, orifice=orifice,
        heave_linear_damping=100, mode_linear_damping=100,
        heave_drag_cd=1.2, heave_drag_area=8,
        pitch_drag_cd=1.2, pitch_drag_area=8,
    )
    result = model.run(
        PMWave(height=1, period=4, phase_file=phase_file),
        dt=0.005, end_time=130, ramp_time=10, radiation_memory=15,
    )
    assert result.case["constraint"]["orifice_force_path"] == "coupled"
    _max_error(result.time, source[:, 0], 1e-10, "OWC time")
    _max_error(result.wave_elevation, _load("wave")[:, 1], 1e-11,
               "OWC incident elevation")
    unwrapped_pitch = np.r_[
        0, np.cumsum((source[1:, 11] + source[:-1, 11]) * 0.005 / 2),
    ]
    early = result.time <= 6
    source_orifice = _load("orifice")
    assert not np.any(source_orifice[source_orifice[:, 0] <= 6, 2])
    assert np.any(source_orifice[source_orifice[:, 0] > 10, 2])
    _max_error(result.bodies["OWC"].position[early, 0],
               source[early, 1], 0.012, "OWC early surge")
    _max_error(result.bodies["OWC"].position[early, 2],
               source[early, 3], 0.025, "OWC early heave")
    _max_error(result.coordinates["pitch_unwrapped"].position[early],
               unwrapped_pitch[early], 0.010, "OWC early physical pitch")
    _max_error(result.flexible_modes["OWC"].position[early, 0],
               flexible[early, 0], 0.0075, "OWC early flexible mode")
    _max_error(result.bodies["OWC"].position[:, 4], source[:, 5],
               0.013, "OWC Euler pitch branch")
    assert np.count_nonzero(
        np.abs(result.bodies["OWC"].position[:, 3] - source[:, 4]) > 1
    ) <= 2
    assert np.count_nonzero(
        np.abs(result.bodies["OWC"].position[:, 5] - source[:, 6]) > 1
    ) <= 2
    pto = result.ptos["orifice"]
    _max_error(pto.absorbed_power, -pto.force * pto.velocity,
               1e-8, "coupled orifice passivity")
    np.testing.assert_array_equal(
        dict(result.raw.extra_outputs)["orifice_compressibility_flag"],
        orifice.evaluate(pto.velocity).compressibility_flag,
    )
    assert np.isfinite(result.bodies["OWC"].position).all()
    assert np.isfinite(result.flexible_modes["OWC"].position).all()


@pytest.mark.skipif(not APPLICATIONS,
                    reason="BEMIO-generated OWC HDF5 absent")
def test_owc_published_force_path_motion(tmp_path):
    """Pair opt-in source force routing, body motion, and PTO over 130 s."""
    hydro = (Path(APPLICATIONS)
             / "OWC/OrificeModel/hydroData/test17a_clean.h5")
    source = _load("OrificeModel_body1")
    flexible = loadmat(
        Path(REFERENCE) / "OWC_ORIFICE_Flex_out.mat",
        simplify_cells=True,
    )["Flex_out"]["signals"]["values"]
    phase_file = tmp_path / "owc_phase.csv"
    np.savetxt(phase_file, _load("components")[:, 3], delimiter=",")
    model = WEC("OWC source force path")
    body = model.body("OWC", hydro, inertia=(99.28, 11.04, 99.2))
    model.floating_gbm(
        body, orifice=OrificePTO(*_load("parameters").ravel()),
        orifice_force_path="published_owc",
        heave_linear_damping=100, mode_linear_damping=100,
        heave_drag_cd=1.2, heave_drag_area=8,
        pitch_drag_cd=1.2, pitch_drag_area=8,
    )
    result = model.run(
        PMWave(height=1, period=4, phase_file=phase_file),
        dt=0.005, end_time=130, ramp_time=10, radiation_memory=15,
    )
    assert result.case["constraint"]["orifice_force_path"] == "published_owc"
    _max_error(result.bodies["OWC"].position[:, 0],
               source[:, 1], 0.175, "source-path rigid surge")
    unwrapped_pitch = np.r_[
        0, np.cumsum((source[1:, 11] + source[:-1, 11]) * 0.005 / 2),
    ]
    _max_error(result.coordinates["pitch_unwrapped"].position,
               unwrapped_pitch, 0.012, "source-path physical pitch")
    _max_error(result.bodies["OWC"].position[:, 2],
               source[:, 3], 0.038, "source-path rigid heave")
    _max_error(result.flexible_modes["OWC"].position[:, 0],
               flexible[:, 0], 0.006, "source-path flexible displacement")
    _max_error(result.bodies["OWC"].velocity[:, 0],
               source[:, 7], 0.008, "source-path rigid surge speed")
    _max_error(result.bodies["OWC"].velocity[:, 2],
               source[:, 9], 0.080, "source-path rigid heave speed")
    _max_error(result.bodies["OWC"].velocity[:, 4],
               source[:, 11], 0.0045, "source-path rigid pitch speed")
    _max_error(result.flexible_modes["OWC"].velocity[:, 0],
               flexible[:, 1], 0.051, "source-path flexible speed")
    source_orifice = _load("orifice")
    source_force = np.interp(result.time, source_orifice[:, 0],
                             source_orifice[:, 1])
    source_power = 1000 * np.interp(result.time, source_orifice[:, 0],
                                    source_orifice[:, 5])
    pto = result.ptos["orifice"]
    _max_error(pto.force, source_force,
               0.06 * np.max(np.abs(source_force)),
               "source-path PTO reaction")
    _max_error(pto.absorbed_power, source_power,
               0.08 * np.max(np.abs(source_power)),
               "source-path PTO absorbed power")
    source_work = np.trapezoid(source_power, result.time)
    model_work = np.trapezoid(pto.absorbed_power, result.time)
    assert source_work > 0 and abs(model_work - source_work) < 0.06 * source_work
    first_ten = result.time <= 10
    _max_error(result.bodies["OWC"].position[first_ten, 2],
               source[first_ten, 3], 0.002, "source-path first-ten heave")
    _max_error(result.flexible_modes["OWC"].position[first_ten, 0],
               flexible[first_ten, 0], 0.0005, "source-path first-ten flex")

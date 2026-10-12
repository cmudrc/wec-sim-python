"""Mapped-coordinate case adapter."""

import numpy as np

from wecsim.components.directLinearGenerator import DirectLinearGenerator
from wecsim.dynamics.case_inputs import (
    _body_number,
    _current_settings,
    _irregular_components_from_case,
    _morison_elements,
    _number,
    _sampled_elevation,
    _section,
)
from wecsim.dynamics.case_response import CaseResponse
from wecsim.dynamics.generalDynamics import BodyMotion, DynamicBody, GeneralizedDynamics
from wecsim.dynamics.hingePitch import hinged_pitch_motion
from wecsim.dynamics.linearCoordinates import build_coordinate_maps, initial_coordinate
from wecsim.dynamics.ptoConnections import build_linear_ptos
from wecsim.environment.irregularWave import (
    imported_spectrum_components,
    synthesize_irregular_response,
)
from wecsim.hydrodynamics.bodyClass import BodyClass
from wecsim.hydrodynamics.morison import (
    no_wave_heave_morison_terms,
    regular_wave_axial_morison_terms,
    regular_wave_heave_morison_terms,
    regular_wave_normal_tangential_morison_terms,
)
from wecsim.hydrodynamics.nonlinearHydro import HeaveMeshHydro
from wecsim.hydrodynamics.passiveYaw import (
    HeldPassiveYawExcitation,
    NearestHeadingExcitation,
    NearestSampledHeadingExcitation,
    PassiveYawExcitation,
    SampledPassiveYawExcitation,
)


def _run_linear_subspace(case, sim, wave, constraint, bodies, hydro,
                         b2b, dt, end_time, ramp_time, rho, g, base_dir):
    """Run mapped coordinates with supported regular, sampled, or no waves."""
    if set(constraint) - {"kind", "initial_coordinate", "initial_speed",
                           "coordinates"}:
        raise ValueError("linear_subspace uses coordinate maps, not joint locations")
    if wave["type"] not in ("regular", "regularCIC", "pm", "jonswap",
                            "spectrumImport", "elevationImport", "none"):
        raise ValueError("linear_subspace supports regular, regularCIC, PM, JONSWAP, spectrumImport, elevationImport, or no waves")
    if b2b and len(set(hydro)) != 1:
        raise ValueError("body-to-body hydrodynamics need one shared HDF5 file")
    if wave["type"] == "none" and set(wave) != {"type"}:
        raise ValueError("no-wave cases have no wave height or period")
    if (wave["type"] in ("pm", "jonswap", "spectrumImport", "elevationImport")
            and any(body.get("mean_drift", "none") != "none" for body in bodies)):
        raise ValueError("irregular linear-subspace mean-drift forcing is not supported")
    components = None
    auxiliary_files = []
    if wave["type"] in ("regular", "regularCIC"):
        allowed = {"type", "height", "period", "direction"}
        if wave["type"] == "regular":
            allowed.add("current")
        if set(wave) - allowed:
            raise ValueError("regular waves use height, period, and direction")
        height = _number(wave.get("height"), "wave.height", nonnegative=True)
        period = _number(wave.get("period"), "wave.period", positive=True)
        direction = _number(wave.get("direction", 0), "wave.direction")
        frequency = 2 * np.pi / period
        if wave["type"] == "regular" and "radiation_memory" in sim:
            raise ValueError("regular-wave linear dynamics use constant radiation")
    elif wave["type"] in ("pm", "jonswap"):
        height = _number(wave.get("height"), "wave.height", positive=True)
        period = _number(wave.get("period"), "wave.period", positive=True)
        components, phase_files = _irregular_components_from_case(
            wave, hydro[0], base_dir,
        )
        auxiliary_files.extend(phase_files)
    elif wave["type"] == "spectrumImport":
        if set(wave) != {"type", "file"} or not isinstance(wave["file"], str) or not wave["file"]:
            raise ValueError("spectrumImport needs one MAT spectrum file")
        spectrum_file = (base_dir / wave["file"]).expanduser().resolve(strict=True)
        components = imported_spectrum_components(hydro[0], spectrum_file)
        auxiliary_files.append(spectrum_file)
    elif wave["type"] == "elevationImport":
        pass  # The sampled record is loaded after the simulation grid is built.
    else:
        if "ramp_time" in sim:
            raise ValueError("ramp_time is inapplicable to no-wave dynamics")
    if wave["type"] != "regular":
        # Radiation convolution samples the HDF5 impulse response on the
        # solver grid; its window must end exactly on a time step.
        memory_time = _number(
            sim.get("radiation_memory", 15), "simulation.radiation_memory",
            positive=True,
        )
        memory_steps = round(memory_time / dt)
        if not np.isclose(memory_steps * dt, memory_time, atol=1e-10):
            raise ValueError("radiation_memory must be a multiple of dt")
        convolution_time = np.arange(memory_steps + 1) * dt
    time = np.arange(round(end_time / dt) + 1) * dt
    if not np.isclose(time[-1], end_time, atol=1e-10):
        raise ValueError("end_time must be an integer multiple of dt")
    imported_elevation = None
    imported_ramp = None
    if wave["type"] == "elevationImport":
        imported_elevation, imported_ramp, wave_path = _sampled_elevation(
            wave, base_dir, time, ramp_time,
        )
        auxiliary_files.append(wave_path)
    body_names = []
    for index, body_spec in enumerate(bodies, start=1):
        _body_number(body_spec, index)
        if "pitch_inertia" in body_spec:
            raise ValueError("linear_subspace uses body.inertia, not pitch_inertia")
        name = body_spec.get("name", f"body{index}")
        if not isinstance(name, str) or not name or name in body_names:
            raise ValueError("body names must be nonempty and unique")
        body_names.append(name)
    loaded_bodies = []
    centers = []
    for index, path in enumerate(hydro, start=1):
        body = BodyClass(str(path))
        body.bodyNumber = index
        body.bodyTotal = len(bodies)
        drift_option = bodies[index - 1].get("mean_drift", "none")
        drift_flags = {"none": 0, "control_surface": 1,
                       "momentum_conservation": 2}
        if not isinstance(drift_option, str) or drift_option not in drift_flags:
            raise ValueError(
                "body.mean_drift must be none, control_surface, or momentum_conservation"
            )
        if drift_option != "none" and wave["type"] not in ("regular", "regularCIC"):
            raise ValueError("mean drift requires regular incident waves")
        body.meanDriftForce = drift_flags[drift_option]
        body.readH5file()
        if drift_option != "none":
            drift_data = body.hydroData["hydro_coeffs"]["mean_drift"]
            if (drift_data.ndim != 3 or drift_data.shape[0] != 6
                    or not np.isfinite(drift_data).all()):
                raise ValueError(
                    f"body{index} HDF5 lacks finite {drift_option} mean-drift coefficients"
                )
        if int(np.asarray(body.dof).item()) != 6:
            raise ValueError("linear_subspace needs six-DOF hydrodynamic bodies")
        center = np.asarray(body.cg, dtype=float).ravel()
        if center.shape != (3,):
            raise ValueError("body center must have three coordinates")
        loaded_bodies.append(body)
        centers.append(center)
    maps, coordinate_names = build_coordinate_maps(
        constraint, bodies, body_names, centers,
    )
    n = maps[0].shape[1]
    exact_pitch_motion = None
    # A pitch coordinate around an offset world pivot also translates the
    # center of gravity. Use exact hinge kinematics for this one-coordinate
    # case instead of treating its map as constant.
    if len(bodies) == n == 1 and "coordinates" in constraint:
        definitions = constraint["coordinates"][0]["motions"]
        if len(definitions) == 1:
            definition = definitions[0]
            pivot = definition.get("pivot", {})
            if (definition["dof"] == "pitch"
                    and float(definition.get("scale", 1)) == 1
                    and "world" in pivot
                    and np.allclose(centers[0][:2], pivot["world"][:2],
                                    rtol=0, atol=1e-12)):
                if any(pto.get("kind") == "linear_actuator"
                       for pto in case.get("ptos", [])):
                    raise ValueError(
                        "exact pitch-pivot motion currently supports rotational PTOs only"
                    )
                exact_pitch_motion = hinged_pitch_motion(
                    centers[0][2] - pivot["world"][2],
                )
    moving_morison = tuple(
        index for index, spec in enumerate(bodies)
        if spec.get("morison_elements")
    )
    moving_elements = ()
    moving_terms = None
    if "current" in wave and not moving_morison:
        raise ValueError("regular-wave current needs a moving Morison element")
    if moving_morison:
        # The verified moving-element force laws use either pure heave or a
        # fixed surge/heave/pitch basis; arbitrary coordinate maps are not
        # yet a valid projection for these Morison terms.
        heave_map = np.zeros((6, 1))
        heave_map[2, 0] = 1
        three_dof_map = np.zeros((6, 3))
        three_dof_map[0, 0] = 1
        three_dof_map[2, 1] = 1
        three_dof_map[4, 2] = 1
        heave_only = n == 1 and np.array_equal(maps[0], heave_map)
        surge_heave_pitch = (n == 3 and np.array_equal(maps[0], three_dof_map)
                             and wave["type"] == "regular")
        if (len(bodies) != 1 or moving_morison != (0,)
                or not (heave_only or surge_heave_pitch)
                or wave["type"] not in ("none", "regular") or b2b
                or bodies[0].get("nonlinear_hydro") is not None
                or bodies[0].get("passive_yaw", False)
                or bodies[0].get("mean_drift", "none") != "none"):
            raise ValueError(
                "moving Morison needs one hydrodynamic body in pure heave "
                "or regular-wave surge/heave/pitch"
            )
        if wave["type"] == "regular" and direction != 0:
            raise ValueError("moving Morison needs zero-heading regular waves")
        if (wave["type"] == "regular"
                and not np.allclose(centers[0][:2], 0, rtol=0, atol=1e-10)):
            raise ValueError("moving Morison needs a body centered at x=y=0")
        if "current" in wave and not surge_heave_pitch:
            raise ValueError("moving Morison current needs regular-wave surge/heave/pitch")
        current_settings = _current_settings(wave)
        moving_elements = _morison_elements(bodies[0]["morison_elements"])
        normal_tangential = [element.axis is not None
                             for element in moving_elements]
        if any(normal_tangential) and not all(normal_tangential):
            raise ValueError("moving Morison elements must use one coefficient mode")
        if any(normal_tangential) and wave["type"] != "regular":
            raise ValueError("normal/tangential moving Morison needs a regular wave")
        if wave["type"] == "regular":
            raw_depth = np.asarray(loaded_bodies[0].hydroData[
                "simulation_parameters"]["water_depth"]).item()
            depth = (np.inf if str(raw_depth).lower() == "infinite"
                     else float(raw_depth))

            if surge_heave_pitch or all(normal_tangential):
                def moving_terms(at_time, coordinate, speed):
                    pose = np.r_[centers[0], np.zeros(3)] + maps[0] @ coordinate
                    velocity = maps[0] @ speed
                    terms = (regular_wave_normal_tangential_morison_terms
                             if all(normal_tangential)
                             else regular_wave_axial_morison_terms)
                    return terms(
                        moving_elements, position=pose, velocity=velocity,
                        time=at_time, wave_height=height, wave_period=period,
                        ramp_time=ramp_time, water_depth=depth, rho=rho, g=g,
                        current_speed=current_settings[0],
                        current_direction=current_settings[1],
                        current_profile=current_settings[2],
                        current_depth=current_settings[3],
                    )
            else:
                def moving_terms(at_time, coordinate, speed):
                    applied, added = regular_wave_heave_morison_terms(
                        moving_elements, center_z=centers[0][2],
                        heave=coordinate[0], speed=speed[0], time=at_time,
                        wave_height=height, wave_period=period,
                        ramp_time=ramp_time, water_depth=depth, rho=rho, g=g,
                    )
                    force = np.zeros(6)
                    force[2] = applied
                    matrix = np.zeros((6, 6))
                    matrix[2, 2] = added
                    return force, matrix
        else:
            def moving_terms(at_time, coordinate, speed):
                applied, added = no_wave_heave_morison_terms(
                    moving_elements, center_z=centers[0][2],
                    heave=coordinate[0], speed=speed[0], rho=rho,
                )
                force = np.zeros(6)
                force[2] = applied
                matrix = np.zeros((6, 6))
                matrix[2, 2] = added
                return force, matrix
        moving_terms(0, np.zeros(n), np.zeros(n))
    passive_indices = [
        index for index, spec in enumerate(bodies)
        if spec.get("passive_yaw", False) is True
    ]
    if any(not isinstance(spec.get("passive_yaw", False), bool)
           for spec in bodies):
        raise ValueError("body.passive_yaw must be a boolean")
    yaw_thresholds = [
        _number(spec.get("passive_yaw_threshold", 0),
                "body.passive_yaw_threshold", nonnegative=True)
        for spec in bodies
    ]
    if any(threshold and (not bodies[index].get("passive_yaw", False)
                          or wave["type"] != "pm")
           for index, threshold in enumerate(yaw_thresholds)):
        raise ValueError("positive passive_yaw_threshold needs PM passive yaw")
    yaw_banks = [spec.get("yaw_heading_bank") for spec in bodies]
    if any(bank is not None and (not bodies[index].get("passive_yaw", False)
                                 or yaw_thresholds[index]
                                 or wave["type"] not in ("regular", "pm"))
           for index, bank in enumerate(yaw_banks)):
        raise ValueError("yaw_heading_bank needs regular or PM passive yaw without a threshold")
    if passive_indices:
        yaw_map = np.zeros((6, 1))
        yaw_map[5, 0] = 1
        if (len(passive_indices) != 1 or n != 1 or b2b
                or wave["type"] not in ("regular", "pm")
                or not np.array_equal(maps[passive_indices[0]], yaw_map)
                or any(np.any(mapping) for index, mapping in enumerate(maps)
                       if index != passive_indices[0])
                or bodies[passive_indices[0]].get("mean_drift", "none") != "none"):
            raise ValueError(
                "passive yaw currently needs one pure-yaw body, stationary others, "
                "regular or PM waves, and independent radiation"
            )
    nonlinear_indices = [
        index for index, spec in enumerate(bodies)
        if spec.get("nonlinear_hydro") is not None
    ]
    if nonlinear_indices:
        heave_map = np.zeros((6, 1))
        heave_map[2, 0] = 1
        if (len(bodies) != 1 or nonlinear_indices != [0] or n != 1
                or not np.array_equal(maps[0], heave_map)
                or wave["type"] not in ("regular", "regularCIC")
                or direction != 0 or b2b
                or bodies[0]["nonlinear_hydro"] != "instantaneous"
                or bodies[0].get("mean_drift", "none") != "none"
                or bodies[0].get("passive_yaw", False)):
            raise ValueError(
                "instantaneous nonlinear hydro currently needs one pure-heave "
                "body and a zero-direction regular or regularCIC wave"
            )
    for spec in bodies:
        if spec.get("nonlinear_hydro") is None and any(
            key in spec for key in ("geometry_file", "drag_coefficient", "drag_area")
        ):
            raise ValueError("mesh geometry and drag settings require nonlinear_hydro")
    initial_q = initial_coordinate(
        constraint.get("initial_coordinate", [0] * n), coordinate_names,
        "constraint.initial_coordinate",
    )
    initial_v = initial_coordinate(
        constraint.get("initial_speed", [0] * n), coordinate_names,
        "constraint.initial_speed",
    )
    if "pto" in case and "ptos" in case:
        raise ValueError("use either a PTO matrix or PTO connections")
    if "pto" in case:
        pto = _section(case["pto"], "pto", {"kind"},
                       {"kind", "stiffness_matrix", "damping_matrix",
                        "equilibrium_coordinate"})
        if pto["kind"] != "linear":
            raise ValueError("linear_subspace requires a linear PTO matrix")
        stiffness = np.asarray(pto.get("stiffness_matrix", np.zeros((n, n))), dtype=float)
        damping = np.asarray(pto.get("damping_matrix", np.zeros((n, n))), dtype=float)
        if (stiffness.shape != (n, n) or damping.shape != (n, n)
                or not np.isfinite(stiffness).all() or not np.isfinite(damping).all()):
            raise ValueError("PTO matrices must be finite N-by-N matrices")
        equilibrium = np.asarray(pto.get("equilibrium_coordinate", [0] * n),
                                 dtype=float)
        if equilibrium.shape != (n,) or not np.isfinite(equilibrium).all():
            raise ValueError("pto.equilibrium_coordinate must have N finite values")
        if np.any(equilibrium) and not np.any(stiffness @ equilibrium):
            raise ValueError("PTO equilibrium offset must produce a spring force")
    else:
        stiffness = np.zeros((n, n))
        damping = np.zeros((n, n))
        equilibrium = np.zeros(n)
    pto_bias = np.zeros(n)
    connections = ()

    dynamic_bodies = []
    nonlinear_models = []
    passive_model = None
    pm_elevation = None
    for index, (body_spec, body, mapping) in enumerate(
            zip(bodies, loaded_bodies, maps), start=1):
        mass_setting = body_spec.get("mass", "equilibrium")
        if mass_setting != "equilibrium":
            mass_setting = _number(mass_setting, "body.mass", positive=True)
        mesh_model = None
        if index - 1 in nonlinear_indices:
            geometry_file = body_spec.get("geometry_file")
            if not isinstance(geometry_file, str) or not geometry_file:
                raise ValueError("instantaneous nonlinear hydro needs geometry_file")
            geometry_path = (base_dir / geometry_file).expanduser().resolve(strict=True)
            cd = _number(body_spec.get("drag_coefficient", 0),
                         "body.drag_coefficient", nonnegative=True)
            drag_area = _number(body_spec.get("drag_area", 0),
                                "body.drag_area", nonnegative=True)
            depth = np.asarray(
                body.hydroData["simulation_parameters"]["water_depth"]
            ).ravel()
            if depth.size != 1 or not np.isfinite(depth[0]) or depth[0] <= 0:
                raise ValueError("nonlinear hydro needs positive HDF5 water depth")
            if not np.allclose(centers[index - 1][:2], 0, rtol=0, atol=1e-10):
                raise ValueError("heave mesh hydro currently needs CG at x=y=0")
            mesh_model = HeaveMeshHydro.from_stl(
                geometry_path, center_z=centers[index - 1][2], rho=rho,
                gravity=g, depth=float(depth[0]), period=period,
                height=height, ramp_time=ramp_time,
                mass=None if mass_setting == "equilibrium" else mass_setting,
                drag_coefficient=cd, drag_area=drag_area,
            )
            if mass_setting == "equilibrium":
                mass_setting = mesh_model.mass
            auxiliary_files.append(geometry_path)
        nonlinear_models.append(mesh_model)
        body.mass = mass_setting
        inertia = np.asarray(body_spec.get("inertia", [0, 0, 0]), dtype=float)
        if (inertia.shape != (3,) or not np.isfinite(inertia).all()
                or np.any(inertia < 0)):
            raise ValueError("body.inertia must have three nonnegative entries")
        body.hydroStiffness = np.zeros((6, 6))
        body.viscDrag = {
            "Drag": np.zeros((6, 6)), "cd": np.zeros(6),
            "characteristicArea": np.zeros(6),
        }
        body.linearDamping = np.zeros((6, 6))
        wave_amp = np.vstack((
            time, imported_elevation if imported_elevation is not None
            else np.zeros_like(time),
        ))
        # BodyClass preprocesses the source HDF5 coefficients for this wave
        # mode. The resulting force blocks feed the common coordinate solver.
        if wave["type"] == "regular":
            body.hydroForcePre(
                frequency, [direction], 1, np.array([0.0]), [], dt, rho, g,
                "regular", wave_amp, index, len(bodies), 0, 0, int(b2b),
            )
        else:
            irf_time = body.hydroData["hydro_coeffs"]["radiation_damping"][
                "impulse_response_fun"]["t"]
            if memory_time > np.max(irf_time) + 1e-10:
                raise ValueError("radiation_memory exceeds the HDF5 kernel")
            regular_memory = wave["type"] == "regularCIC"
            irregular = wave["type"] in ("pm", "jonswap", "spectrumImport")
            body.hydroForcePre(
                (frequency if regular_memory else
                 components.omega if irregular else []),
                ([direction] if regular_memory else
                 components.directions if irregular else [0]),
                len(convolution_time), convolution_time,
                len(components.omega) if irregular else [],
                dt, rho, g, ("regularCIC" if regular_memory else
                             "spectrumImport" if wave["type"] == "spectrumImport" else
                             "elevationImport" if wave["type"] == "elevationImport" else
                             "irregular" if irregular else "noWaveCIC"),
                wave_amp,
                index, len(bodies), 0, 0, int(b2b),
            )
        physical_mass = float(np.asarray(body.mass).item())
        rigid_mass = np.diag([physical_mass] * 3 + inertia.tolist())
        hydro_force = body.hydroForce
        # A body-to-body HDF5 contains cross-body added-mass/damping blocks.
        # Independent bodies instead contribute only their diagonal block.
        if b2b:
            added_mass = tuple(
                np.asarray(hydro_force["fAddedMass"])[:, 6*j:6*(j+1)]
                for j in range(len(bodies))
            )
            if wave["type"] == "regular":
                radiation_damping = tuple(
                    np.asarray(hydro_force["fDamping"])[:, 6*j:6*(j+1)]
                    for j in range(len(bodies))
                )
        else:
            added_mass = [np.zeros((6, 6)) for _ in bodies]
            added_mass[index - 1] = np.asarray(hydro_force["fAddedMass"])
            if wave["type"] == "regular":
                radiation_damping = [np.zeros((6, 6)) for _ in bodies]
                radiation_damping[index - 1] = np.asarray(hydro_force["fDamping"])
        if wave["type"] != "regular":
            radiation_damping = tuple(np.zeros((6, 6)) for _ in bodies)
            kernel = np.asarray(hydro_force["irkb"])
            if not b2b and len(bodies) > 1:
                independent_kernel = np.zeros((len(kernel), 6, 6 * len(bodies)))
                independent_kernel[:, :, 6 * (index - 1):6 * index] = kernel
                kernel = independent_kernel
        else:
            kernel = None
        center = centers[index - 1]

        if exact_pitch_motion is not None:
            motion = exact_pitch_motion
        else:
            def motion(q, v, *, mapping=mapping):
                return BodyMotion(mapping @ q, mapping, np.zeros(6))

        if wave["type"] in ("regular", "regularCIC"):
            real = np.asarray(hydro_force["fExt"]["re"])
            imaginary = np.asarray(hydro_force["fExt"]["im"])
            drift = np.asarray(hydro_force["fExt"]["md"])

            def excitation(at_time, *, real=real, imaginary=imaginary,
                           drift=drift):
                ramp = (1.0 if ramp_time == 0 or at_time >= ramp_time
                        else (1 - np.cos(np.pi * at_time / ramp_time)) / 2)
                # HDF5 complex excitation follows re*cos(omega*t)
                # minus im*sin(omega*t); mean drift scales with amplitude squared.
                return (height / 2 * ramp * (
                    real * np.cos(frequency * at_time)
                    - imaginary * np.sin(frequency * at_time)
                ) + (height / 2)**2 * ramp * drift)
        elif wave["type"] == "elevationImport":
            sampled_force = np.asarray(hydro_force["userDefinedFe"])
            if wave.get("reapply_force_ramp", False):
                sampled_force = sampled_force * imported_ramp[:, None]

            def excitation(at_time, *, sampled_force=sampled_force):
                if len(time) == 1:
                    return sampled_force[0]
                sample = min(max(at_time / dt, 0.0), float(len(time) - 1))
                left = min(int(sample), len(time) - 2)
                fraction = sample - left
                return ((1 - fraction) * sampled_force[left]
                        + fraction * sampled_force[left + 1])
        elif wave["type"] in ("pm", "jonswap", "spectrumImport"):
            if body_spec.get("passive_yaw", False):
                def excitation(at_time):
                    return np.zeros(6)
            else:
                incident = synthesize_irregular_response(
                    hydro[index - 1], components, dt=dt, end_time=end_time,
                    ramp_time=ramp_time, body_number=index, rho=rho, g=g,
                    excitation_interpolation=wave.get(
                        "excitation_interpolation", "linear",
                    ),
                )
                if pm_elevation is None:
                    pm_elevation = incident.elevation
                sampled_force = incident.excitation_force

                def excitation(at_time, *, sampled_force=sampled_force):
                    sample = round(at_time / dt)
                    return sampled_force[sample]
        else:
            def excitation(at_time):
                return np.zeros(6)

        state_excitation = None
        if body_spec.get("passive_yaw", False):
            if wave["type"] == "pm":
                passive_model = SampledPassiveYawExcitation.from_hydro_data(
                    body.hydroData, components, dt=dt, end_time=end_time,
                    ramp_time=ramp_time, rho=rho, g=g,
                )
                pm_elevation = passive_model.elevation
                if yaw_banks[index - 1] is not None:
                    passive_model = NearestSampledHeadingExcitation(
                        passive_model, yaw_banks[index - 1],
                    )
                if yaw_thresholds[index - 1]:
                    passive_model = HeldPassiveYawExcitation(
                        passive_model, yaw_thresholds[index - 1],
                    )
            else:
                passive_model = PassiveYawExcitation.from_hydro_data(
                    body.hydroData, omega=frequency,
                    incident_direction=direction, amplitude=height / 2,
                    ramp_time=ramp_time, rho=rho, g=g,
                    spline_frequency=yaw_banks[index - 1] is not None,
                )
                if yaw_banks[index - 1] is not None:
                    passive_model = NearestHeadingExcitation(
                        passive_model, yaw_banks[index - 1],
                    )

            if isinstance(passive_model, (HeldPassiveYawExcitation,
                                          NearestSampledHeadingExcitation)):
                state_excitation = passive_model
            else:
                def state_excitation(at_time, coordinate, speed, *,
                                     model=passive_model):
                    return model.force(at_time, coordinate[0])
        if mesh_model is not None:
            def state_excitation(at_time, coordinate, speed, *,
                                 model=mesh_model, linear=excitation):
                buoyancy, fk, drag = model.forces(
                    at_time, coordinate[0], speed[0]
                )
                result = linear(at_time).copy()
                result[2] += buoyancy + fk + drag
                return result

        state_inertia = None
        if moving_elements:
            def state_inertia(at_time, coordinate, speed, *, terms=moving_terms):
                return terms(at_time, coordinate, speed)

        dynamic_bodies.append(DynamicBody(
            rigid_mass=rigid_mass,
            added_mass=tuple(added_mass),
            damping=tuple(radiation_damping),
            restoring=(np.zeros((6, 6)) if mesh_model is not None else
                       np.asarray(hydro_force["linearHydroRestCoef"])),
            static_force=(np.zeros(6) if mesh_model is not None else
                          np.array([
                              0, 0, (rho * float(np.asarray(body.dispVol).item())
                                     - physical_mass) * g, 0, 0, 0,
                          ])),
            reference_position=np.r_[center, np.zeros(3)],
            motion=motion,
            excitation=excitation,
            radiation_kernel=kernel,
            state_excitation=state_excitation,
            state_inertia=state_inertia,
        ))
    if "ptos" in case:
        connections, stiffness, damping, pto_bias = build_linear_ptos(
            case["ptos"], maps, centers, body_names, coordinate_names,
            allow_linear_generator=True,
        )
    generator_specs = {
        spec["name"]: spec["linear_generator"]
        for spec in case.get("ptos", []) if "linear_generator" in spec
    }
    if generator_specs and wave["type"] != "regular":
        raise ValueError("linear generator currently needs regular-wave constant radiation")
    linear_generators = []
    for connection in connections:
        if connection.name not in generator_specs:
            continue
        name = f"pto.{connection.name}.linear_generator"
        spec = _section(
            generator_specs[connection.name], name,
            {"stator_resistance", "friction", "pole_pitch", "magnet_flux",
             "inductance", "load_resistance"},
            {"stator_resistance", "friction", "pole_pitch", "magnet_flux",
             "inductance", "load_resistance", "initial_angle",
             "initial_flux_d", "initial_flux_q"},
        )
        settings = {
            key: _number(value, f"{name}.{key}")
            for key, value in spec.items() if value is not None
        }
        linear_generators.append((connection, DirectLinearGenerator(**settings)))
    linear_generators = tuple(linear_generators)
    controlled_connections = tuple(
        connection for connection in connections if connection.control is not None
    )
    system = GeneralizedDynamics(
        tuple(dynamic_bodies), n, pto_stiffness=stiffness,
        pto_damping=damping, pto_equilibrium=equilibrium,
        pto_bias=pto_bias,
        controlled_ptos=controlled_connections,
        linear_generators=linear_generators,
    )
    response = system.integrate(
        dt=dt, end_time=end_time,
        initial_coordinate=initial_q, initial_speed=initial_v,
    )
    if wave["type"] in ("regular", "regularCIC"):
        ramp = np.ones(len(response.time))
        if ramp_time > 0:
            early = response.time < ramp_time
            ramp[early] = (1 - np.cos(np.pi * response.time[early] / ramp_time)) / 2
        elevation = height / 2 * ramp * np.cos(frequency * response.time)
    else:
        elevation = (imported_elevation if imported_elevation is not None
                     else pm_elevation)
    coordinate_outputs = tuple(
        output for index, name in enumerate(coordinate_names)
        for output in (
            (f"coordinate_{name}_position", response.coordinate[:, index]),
            (f"coordinate_{name}_velocity", response.speed[:, index]),
        )
    ) if "coordinates" in constraint else ()
    controlled_forces = (
        {connection.name: response.controlled_pto_force[:, index]
         for index, connection in enumerate(controlled_connections)}
        if controlled_connections else {}
    )
    generator_forces = (
        {connection.name: response.linear_generator_force[:, index]
         for index, (connection, _) in enumerate(linear_generators)}
        if linear_generators else {}
    )
    pto_outputs = tuple(
        output for connection in connections
        for output in connection.outputs(
            response.coordinate, response.speed,
            force_override=(controlled_forces | generator_forces).get(connection.name),
        )
    )
    generator_outputs = []
    for index, (connection, generator) in enumerate(linear_generators):
        states = response.linear_generator_state[:, index]
        stroke_speed = response.speed @ connection.stroke_jacobian
        electrical_power = np.empty(len(response.time))
        phase_current = np.empty((len(response.time), 3))
        phase_voltage = np.empty_like(phase_current)
        for step, (speed, state) in enumerate(zip(stroke_speed, states)):
            signals = generator.signals(float(speed), state)
            electrical_power[step] = signals.electrical_power
            phase_current[step] = signals.phase_current
            phase_voltage[step] = signals.phase_voltage
        prefix = f"pto_{connection.name}_generator_"
        generator_outputs.extend((
            (prefix + "flux_d", states[:, 0]),
            (prefix + "flux_q", states[:, 1]),
            (prefix + "angle", states[:, 2]),
            (prefix + "friction_force", generator.friction * stroke_speed),
            (prefix + "electrical_power", electrical_power),
            (prefix + "phase_current", phase_current),
            (prefix + "phase_voltage", phase_voltage),
        ))
    generalized_pto = (-(response.coordinate - equilibrium) @ stiffness.T
                       - response.speed @ damping.T + pto_bias)
    # Controlled and generator forces are history-dependent, so add their
    # sampled scalar forces after the passive matrix contribution.
    for connection in controlled_connections:
        generalized_pto += np.outer(
            controlled_forces[connection.name], connection.stroke_jacobian,
        )
    for connection, _ in linear_generators:
        generalized_pto += np.outer(
            generator_forces[connection.name], connection.stroke_jacobian,
        )
    drift_outputs = tuple(
        output
        for index, (spec, body, dynamic) in enumerate(
            zip(bodies, loaded_bodies, dynamic_bodies), start=1,
        )
        if spec.get("mean_drift", "none") != "none"
        for output in (
            (f"body{index}_mean_drift_force",
             (height / 2)**2 * ramp[:, None]
             * np.asarray(body.hydroForce["fExt"]["md"])[None, :]),
            (f"body{index}_excitation_force",
             np.stack([dynamic.excitation(t) for t in response.time])),
        )
    ) if wave["type"] in ("regular", "regularCIC") else ()
    imported_outputs = tuple(
        (f"body{index}_excitation_force",
         (np.asarray(body.hydroForce["userDefinedFe"])
          * (imported_ramp[:, None] if wave.get("reapply_force_ramp", False)
             else 1)))
        for index, body in enumerate(loaded_bodies, start=1)
    ) if wave["type"] == "elevationImport" else ()
    passive_outputs = ()
    if passive_model is not None:
        body_index = passive_indices[0] + 1
        if isinstance(passive_model, HeldPassiveYawExcitation):
            excitation_history = np.asarray(passive_model.force_history)
        else:
            excitation_history = np.stack([
                passive_model.force(t, angle)
                for t, angle in zip(response.time, response.coordinate[:, 0])
            ])
        passive_outputs = ((
            f"body{body_index}_excitation_force",
            excitation_history,
        ),)
    nonlinear_outputs = ()
    if nonlinear_indices:
        model = nonlinear_models[0]
        components = np.array([
            model.forces(t, q[0], v[0])
            for t, q, v in zip(response.time, response.coordinate, response.speed)
        ])
        nonlinear_outputs = (
            ("body1_buoyancy_minus_weight", components[:, 0]),
            ("body1_nonlinear_fk_correction", components[:, 1]),
            ("body1_quadratic_drag_force", components[:, 2]),
        )
    morison_outputs = ()
    if moving_elements:
        # Added mass entered the acceleration solve as inertia. Reconstruct
        # the reported physical Morison force by removing that acceleration
        # term from the applied-force part.
        force = np.stack([
            applied - added @ (maps[0] @ a)
            for (applied, added), a in zip(
                (moving_terms(t, q, v)
                 for t, q, v in zip(response.time, response.coordinate,
                                    response.speed)),
                response.acceleration,
            )
        ])
        morison_outputs = ((f"morison_force_{body_names[0]}", force),)
    return CaseResponse(
        response.time, response.body_position, response.body_velocity,
        hydro, wave_elevation=elevation,
        auxiliary_files=tuple(auxiliary_files),
        pto_generalized_force=(
            generalized_pto if "pto" in case or "ptos" in case else None
        ),
        extra_outputs=(coordinate_outputs + pto_outputs + tuple(generator_outputs) + drift_outputs
                       + imported_outputs
                       + passive_outputs + nonlinear_outputs + morison_outputs),
    )

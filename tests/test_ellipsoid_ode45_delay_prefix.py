"""Pin the active added-mass feedback in the published ode45 ellipsoid."""

from pathlib import Path

import numpy as np

from wecsim.hydrodynamics.bodyClass import BodyClass
from wecsim.hydrodynamics.nonlinearHydro import HeaveMeshHydro


ROOT = Path(__file__).resolve().parents[1]
FIXTURE = (ROOT / "tests/fixtures/matlab"
           / "ellipsoid_ode45_delay_2s/signals.npz")
CLOCK = (ROOT / "tests/fixtures/matlab/ellipsoid_ode45_nonlinear_clock"
         / "source_maxstep_common_grid.npz")
HYDRO = ROOT / "tests/test_objects/test_bodyclass/testData/hydroData/ellipsoid.h5"
MESH = ROOT / "tests/test_objects/test_bodyclass/testData/body_10_test/elipsoid.stl"


def test_published_ellipsoid_added_mass_uses_delayed_acceleration():
    with np.load(FIXTURE, allow_pickle=False) as source:
        time = source["time"]
        acceleration = source["delay_input_heave"]
        delayed = source["delay_output_heave"]
        body_acceleration = source["body_acceleration_heave"]
        added_force = source["applied_added_mass_heave"]
        total_force = source["total_hydrodynamic_force_heave"]
        pto_force = source["pto_force_heave"]
        adjusted_mass = float(source["source_adjusted_mass"])
        applied_coefficient = float(
            source["source_applied_added_mass_coefficient"]
        )

    np.testing.assert_allclose(time, np.arange(41) * .05,
                               rtol=0, atol=1e-13)
    np.testing.assert_array_equal(acceleration, body_acceleration)
    np.testing.assert_array_equal(delayed[:2], [0., 0.])
    # The 0.1 microsecond Transport Delay extrapolates from the two prior
    # major-step accelerations. The small delay correction matters when
    # checking this source mechanism at its actual precision.
    expected = (2 * acceleration[1:-1] - acceleration[:-2]
                - (1e-7 / .05) * (acceleration[1:-1]
                                   - acceleration[:-2]))
    assert np.max(np.abs(delayed[2:] - expected)) < 2e-16
    assert np.max(np.abs(added_force - applied_coefficient * delayed)) < 1e-8
    assert np.max(np.abs(adjusted_mass * body_acceleration
                         - total_force - pto_force)) < 1e-7


def test_python_force_law_closes_published_ellipsoid_rhs_on_source_states():
    """Rebuild the source force balance from pinned mesh and BEM inputs."""
    with np.load(CLOCK, allow_pickle=False) as source:
        time = source["regular_published_time"][:41]
        heave = source["regular_published_heave"][:41] + 2
        speed = source["regular_published_speed"][:41]
    with np.load(FIXTURE, allow_pickle=False) as source:
        np.testing.assert_allclose(time, source["time"], rtol=0, atol=1e-13)
        delayed = source["delay_output_heave"]
        expected = source["body_acceleration_heave"]
        source_mass = float(source["source_adjusted_mass"])
        source_coefficient = float(source["source_applied_added_mass_coefficient"])

    mesh = HeaveMeshHydro.from_stl(
        MESH, center_z=-2, rho=1025, gravity=9.81, depth=70,
        period=6, height=4, ramp_time=50, mass=None,
        drag_coefficient=1, drag_area=np.pi * 25,
    )
    body = BodyClass(str(HYDRO))
    body.bodyNumber = body.bodyTotal = 1
    body.readH5file()
    body.mass = mesh.mass
    body.inertia = [1.375264e6, 1.375264e6, 1.341721e6]
    body.hydroStiffness = np.zeros((6, 6))
    body.viscDrag = {"Drag": np.zeros((6, 6)), "cd": np.zeros(6),
                     "characteristicArea": np.zeros(6)}
    body.linearDamping = np.zeros((6, 6))
    frequency = 2 * np.pi / 6
    body.hydroForcePre(frequency, [0], 1, np.array([0.]), [], .05,
                       1025, 9.81, "regular", np.zeros((2, 2)), 1, 1, 0, 0, 0)
    hydro = body.hydroForce
    added_mass = np.asarray(hydro["fAddedMass"])
    adjusted_mass = mesh.mass + 2 * np.trace(added_mass[:3, :3])
    # The source moves the difference between this mass and the ordinary
    # heave mass to a delayed feedback force to break its algebraic loop.
    feedback = adjusted_mass - (mesh.mass + added_mass[2, 2])
    np.testing.assert_allclose(adjusted_mass, source_mass, rtol=0, atol=1e-8)
    np.testing.assert_allclose(feedback, -source_coefficient, rtol=0, atol=1e-8)

    real = np.asarray(hydro["fExt"]["re"])[2]
    imaginary = np.asarray(hydro["fExt"]["im"])[2]
    damping = float(np.asarray(hydro["fDamping"])[2, 2]) + 1_200_000
    reconstructed = np.empty_like(expected)
    for index, at_time in enumerate(time):
        prior = max(0, index - 1)
        buoyancy = mesh.forces(time[prior], heave[prior], speed[prior])[0]
        _, froude_krylov, drag = mesh.forces(
            at_time, heave[index], speed[index],
        )
        linear = 2 * mesh.ramp(at_time) * (
            real * np.cos(frequency * at_time)
            - imaginary * np.sin(frequency * at_time)
        )
        force = (linear + buoyancy + froude_krylov + drag
                 - damping * speed[index])
        reconstructed[index] = (force + feedback * delayed[index]) / adjusted_mass
    assert np.max(np.abs(reconstructed - expected)) < 1e-11

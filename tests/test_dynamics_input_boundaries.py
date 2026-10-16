"""Bad physical blocks must fail before the shared solver advances a WEC."""

from dataclasses import replace
from types import SimpleNamespace

import numpy as np
import pytest

from wecsim.dynamics.generalDynamics import (
    BodyMotion, DynamicBody, GeneralizedDynamics,
)


def _body():
    jacobian = np.zeros((6, 1))
    jacobian[2, 0] = 1
    mass = np.zeros((6, 6))
    mass[2, 2] = 2
    zero_matrix = np.zeros((6, 6))

    def motion(coordinate, speed):
        return BodyMotion(
            jacobian[:, 0] * coordinate[0], jacobian, np.zeros(6),
        )

    return DynamicBody(
        mass, (zero_matrix,), (zero_matrix,), zero_matrix,
        np.zeros(6), np.zeros(6), motion, lambda time: np.zeros(6),
    )


@pytest.mark.parametrize("change, error", (
    ({"rigid_mass": np.zeros((5, 5))}, "finite 6x6 matrices"),
    ({"restoring": np.full((6, 6), np.nan)}, "finite 6x6 matrices"),
    ({"added_mass": (np.zeros((6, 5)),)}, "finite 6x6 matrices"),
    ({"damping": ()}, "one hydrodynamic block"),
    ({"static_force": np.zeros(5)}, "finite six-vectors"),
    ({"reference_position": np.full(6, np.inf)}, "finite six-vectors"),
    ({"radiation_kernel": np.zeros((2, 6, 5))}, "radiation kernel"),
    ({"state_excitation": 1}, "state_excitation must be callable"),
    ({"state_inertia": 1}, "state_inertia must be callable"),
))
def test_device_rejects_malformed_body_physics(change, error):
    with pytest.raises((ValueError, TypeError), match=error):
        GeneralizedDynamics((replace(_body(), **change),), 1)


@pytest.mark.parametrize("settings, error", (
    ({"pto_stiffness": np.zeros((2, 2))}, "PTO matrices"),
    ({"pto_damping": np.full((1, 1), np.nan)}, "PTO matrices"),
    ({"pto_equilibrium": np.zeros(2)}, "PTO matrices"),
    ({"radiation_discretization": "unknown"}, "radiation_discretization"),
    ({"radiation_discretization": "fir"}, "FIR radiation needs a kernel"),
    ({"added_mass_delay": .01}, "added-mass delay needs"),
    ({"nonlinear_force": 1}, "nonlinear_force must be callable"),
    ({"controlled_ptos": (SimpleNamespace(control=None),)},
     "controlled PTO needs a finite"),
))
def test_device_rejects_inconsistent_coordinate_force_laws(settings, error):
    with pytest.raises((ValueError, TypeError), match=error):
        GeneralizedDynamics((_body(),), 1, **settings)


@pytest.mark.parametrize("force_kind, force", (
    ("excitation", np.full(6, np.nan)),
    ("nonlinear", np.array([np.nan])),
    ("applied", np.array([np.nan])),
))
def test_device_rejects_nonfinite_runtime_forces(force_kind, force):
    body = _body()
    options = {}
    if force_kind == "excitation":
        body = replace(body, excitation=lambda time: force)
    elif force_kind == "nonlinear":
        options["nonlinear_force"] = lambda coordinate, speed: force
    system = GeneralizedDynamics((body,), 1, **options)
    with pytest.raises(ValueError, match="finite"):
        system.acceleration(
            0, np.zeros(1), np.zeros(1),
            **({"applied_force": force} if force_kind == "applied" else {}),
        )


def test_device_rejects_singular_constrained_mass():
    body = replace(_body(), rigid_mass=np.zeros((6, 6)))
    with pytest.raises(ValueError, match="singular mass matrix"):
        GeneralizedDynamics((body,), 1).acceleration(0, np.zeros(1), np.zeros(1))

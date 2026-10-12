"""Compare the published Sphere reactive controller with paired MATLAB output.

The pinned Simulink controller applies -Kp * heave_velocity - Ki * heave_position
to the translational PTO, with Ki negative in this published example. The
signed matrix in the generalized-coordinate runner represents that force.
"""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim.dynamics.caseDynamics import run_case


SPHERE_H5 = os.environ.get("WEC_SIM_SPHERE_H5")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
pytestmark = pytest.mark.skipif(
    not (SPHERE_H5 and REFERENCE),
    reason="paired MATLAB Sphere reactive-controller output not provided",
)

KP = 49_181.0
KI = -573_350.0


def _max_error(actual, expected, limit, label):
    assert actual.shape == expected.shape, label
    assert np.isfinite(expected).all(), label
    error = np.max(np.abs(actual - expected))
    assert error < limit, f"{label}: max error {error:.6g} exceeds {limit}"


def test_published_sphere_reactive_pi_against_matlab():
    hydro = str(Path(SPHERE_H5).resolve())
    case = {
        "simulation": {"dt": 0.02, "end_time": 200, "ramp_time": 50},
        "wave": {"type": "regular", "height": 2.5, "period": 9.6664},
        "bodies": [{
            "hydro_file": hydro, "name": "sphere", "mass": "equilibrium",
            "inertia": [20_907_301, 21_306_090.66, 37_085_481.11],
        }],
        "constraint": {"kind": "linear_subspace", "coordinates": [{
            "name": "heave", "motions": [{"body": "sphere", "dof": "heave"}],
        }]},
        "pto": {"kind": "linear", "stiffness_matrix": [[KI]],
                "damping_matrix": [[KP]]},
    }
    response = run_case(case)
    reference = Path(REFERENCE)
    body = np.loadtxt(
        reference / "Sphere_Reactive_PI_Reactive_PI_body1.csv", delimiter=",",
    )
    controller = np.loadtxt(
        reference / "Sphere_Reactive_PI_controller.csv", delimiter=",",
    )
    assert body.shape == (10_001, 25)
    assert controller.shape == (10_001, 13)
    np.testing.assert_allclose(response.time, body[:, 0], rtol=0, atol=1e-10)
    np.testing.assert_allclose(response.time, controller[:, 0], rtol=0, atol=1e-10)

    outputs = dict(response.extra_outputs)
    displacement = outputs["coordinate_heave_position"]
    speed = outputs["coordinate_heave_velocity"]
    force = -KP * speed - KI * displacement
    power = force * speed
    _max_error(response.body_position[:, 0, 2], body[:, 3], 0.010,
               "Sphere heave position")
    _max_error(response.body_velocity[:, 0, 2], body[:, 9], 0.006,
               "Sphere heave velocity")
    _max_error(force, controller[:, 3], 5_500, "reactive controller force")
    _max_error(power, controller[:, 9], 40_000, "reactive controller power")

    # The saved MATLAB controller must also follow its documented gains.
    matlab_force = -KP * body[:, 9] - KI * (body[:, 3] - body[0, 3])
    _max_error(matlab_force, controller[:, 3], 1e-4,
               "MATLAB controller gain relation")
    _max_error(controller[:, 3] * body[:, 9], controller[:, 9], 1e-3,
               "MATLAB controller power relation")

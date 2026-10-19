"""Separate WaveStar plant integration from autonomous NMPC feedback.

The sole source-driven runtime input is the published controller's sampled
requested torque. Body motion, axial PTO force, and radiation are advanced by
Python from the pinned hydrodynamic data and incident-wave phases.
"""

from pathlib import Path

import numpy as np

from wecsim import WaveStarNmpcActuator
from wecsim.dynamics.wavestar import run_wavestar_published
from wecsim.environment.irregularWave import jonswap_equal_energy_components


ROOT = Path(__file__).parent / "fixtures/matlab"
HYDRO = ROOT / "wavestar_fixed_input_step/wavestar.h5"
SOURCE = ROOT / "wavestar_nmpc_full_step/baseline.npz"
PYTHON_COMMAND_SOURCE = ROOT / "wavestar_nmpc_full_step/python_command_source.npz"


class _PublishedTorqueReplay:
    """Apply the logged request on its own tick; the actuator supplies its lag."""

    def __init__(self, command, *, plant_dt=.05):
        self.plant_dt = plant_dt
        self.stride = round(.05 / plant_dt)
        if not np.isclose(self.stride * plant_dt, .05, rtol=0, atol=1e-12):
            raise ValueError("plant step must divide the published controller step")
        self.command = command
        self.actuator = WaveStarNmpcActuator()
        self.index = 0

    def step(self, stroke, *, noise=0., dropout=False):
        assert noise == 0. and not dropout
        if self.index % self.stride == 0:
            force = self.actuator.step(
                self.command[self.index // self.stride], stroke,
            )
        else:
            force = self.actuator.force_at(stroke)
        self.index += 1
        return force

    def force_at_stroke(self, stroke):
        return self.actuator.force_at(stroke)


def test_published_225s_plant_with_prescribed_nmpc_requests():
    with np.load(SOURCE, allow_pickle=False) as source:
        body = source["body"]
        pto = source["pto"]
        command = source["command"]
        command_time = source["command_time"]
        phase = source["components"][:, 3, None]
    assert body.shape == (4501, 25) and pto.shape == (4501, 43)
    assert command.shape == command_time.shape == (4501,)
    assert np.isfinite(body).all() and np.isfinite(pto).all()
    assert np.isfinite(command).all() and np.isfinite(phase).all()
    np.testing.assert_allclose(command_time, body[:, 0], rtol=0, atol=1e-10)
    assert np.max(np.abs(command[command_time >= 15])) > 1

    sea = jonswap_equal_energy_components(
        HYDRO, significant_height=.1042, peak_period=1.836,
        directions=np.array([0.]), spreading=np.array([1.]),
        gamma=3.3, phase=phase,
    )
    response = run_wavestar_published(
        HYDRO, sea, dt=.05, end_time=225, ramp_time=25,
        g=9.80665, solver="wecsim_ode8",
        pto_controller=_PublishedTorqueReplay(command),
    )
    np.testing.assert_allclose(response.time, body[:, 0], rtol=0, atol=1e-10)
    assert np.max(np.abs(response.angle - body[:, 5])) < 3e-4
    assert np.max(np.abs(response.angular_speed - body[:, 11])) < 2e-3
    assert np.max(np.abs(response.pto_stroke - pto[:, 3])) < 6e-5
    assert np.max(np.abs(response.pto_force - pto[:, 33])) < .005
    python_work = np.trapezoid(
        response.pto_force * response.pto_speed, response.time,
    )
    source_work = np.trapezoid(pto[:, 33] * pto[:, 9], pto[:, 0])
    assert abs(python_work - source_work) / abs(source_work) < .001


def test_225s_source_plant_on_independent_python_nmpc_requests():
    """Check the reciprocal plant solve on the Python controller's commands.

    MATLAB recomputes its own hydrodynamics and plant trajectory; only the
    active requested motor torque is prescribed. This is a component gate,
    not an autonomous closed-loop trajectory comparison.
    """
    with np.load(PYTHON_COMMAND_SOURCE, allow_pickle=False) as source:
        time = source["time"]
        command = source["python_command"]
        logged_command = source["source_command"]
        pitch = source["pitch"]
        speed = source["pitch_speed"]
        stroke = source["pto_stroke"]
        pto_speed = source["pto_speed"]
        force = source["pto_force"]
    with np.load(SOURCE, allow_pickle=False) as baseline:
        phase = baseline["components"][:, 3, None]

    assert all(values.shape == (4501,) for values in (
        time, command, logged_command, pitch, speed, stroke, pto_speed, force,
    ))
    assert all(np.isfinite(values).all() for values in (
        time, command, logged_command, pitch, speed, stroke, pto_speed, force,
    ))
    active = time >= 15
    assert np.count_nonzero(active) == 4201
    np.testing.assert_array_equal(logged_command[active], command[active])
    sea = jonswap_equal_energy_components(
        HYDRO, significant_height=.1042, peak_period=1.836,
        directions=np.array([0.]), spreading=np.array([1.]),
        gamma=3.3, phase=phase,
    )
    response = run_wavestar_published(
        HYDRO, sea, dt=.05, end_time=225, ramp_time=25,
        g=9.80665, solver="wecsim_ode8",
        pto_controller=_PublishedTorqueReplay(command),
    )
    np.testing.assert_allclose(response.time, time, rtol=0, atol=1e-10)
    assert np.max(np.abs(response.angle - pitch)) < 3e-4
    assert np.max(np.abs(response.angular_speed - speed)) < 2e-3
    assert np.max(np.abs(response.pto_stroke - stroke)) < 6e-5
    # MATLAB's precontrol resistance law differs slightly; the prescribed
    # NMPC torque begins at 15 s, so gate axial force on that interval.
    assert np.max(np.abs(response.pto_force[active] - force[active])) < .005
    python_work = np.trapezoid(
        response.pto_force * response.pto_speed, response.time,
    )
    source_work = np.trapezoid(force * pto_speed, time)
    assert abs(python_work - source_work) / abs(source_work) < .001


def test_physical_plant_gap_persists_with_published_requests():
    """Separate physical plant integration from autonomous NMPC feedback.

    This is a negative published-case diagnostic. The opt-in source-ode8
    plant above pairs the source on the same requests; these two physical
    runs instead retain implicit added mass and refine only Python's step.
    """
    with np.load(SOURCE, allow_pickle=False) as source:
        body = source["body"]
        pto = source["pto"]
        command = source["command"]
        phase = source["components"][:, 3, None]
    sea = jonswap_equal_energy_components(
        HYDRO, significant_height=.1042, peak_period=1.836,
        directions=np.array([0.]), spreading=np.array([1.]),
        gamma=3.3, phase=phase,
    )
    source_work = np.trapezoid(pto[:, 33] * pto[:, 9], pto[:, 0])
    assert abs(source_work) > 100
    responses = []
    for step in (.025, .0125):
        stride = round(.05 / step)
        response = run_wavestar_published(
            HYDRO, sea, dt=step, end_time=225, ramp_time=25,
            g=9.80665, output_stride=stride,
            pto_controller=_PublishedTorqueReplay(command, plant_dt=step),
        )
        np.testing.assert_allclose(response.time, body[:, 0],
                                   rtol=0, atol=1e-10)
        assert np.isfinite(response.angle).all()
        assert np.isfinite(response.pto_force).all()
        # The same requests give nearly the same actuator force, while the
        # physical plant still misses the published source trajectory.
        assert np.max(np.abs(response.pto_force - pto[:, 33])) < 1
        pitch_gap = np.max(np.abs(response.angle - body[:, 5]))
        work = np.trapezoid(response.pto_force * response.pto_speed,
                            response.time)
        work_gap = abs(work - source_work) / abs(source_work)
        assert .03 < pitch_gap < .05
        assert .05 < work_gap < .10
        responses.append((response, work))

    coarse, coarse_work = responses[0]
    fine, fine_work = responses[1]
    # A stable Python refinement alongside a large source gap rules out
    # Python's plant step and the optimizer as explanations for this result.
    assert np.max(np.abs(coarse.angle - fine.angle)) < 1e-4
    assert np.max(np.abs(coarse.angular_speed - fine.angular_speed)) < 1e-3
    assert abs(coarse_work - fine_work) / abs(source_work) < .001

"""Isolate the pinned ode4 and radiation extrapolation path for held yaw."""

import os
from pathlib import Path

import numpy as np
import pytest
from scipy.interpolate import CubicSpline, PchipInterpolator

from wecsim.environment.irregularWave import pm_equal_energy_components
from wecsim.hydrodynamics.bodyClass import BodyClass
from wecsim.hydrodynamics.passiveYaw import (
    HeldPassiveYawExcitation, SampledPassiveYawExcitation,
)


APPLICATIONS = os.environ.get("WEC_SIM_APPLICATIONS_DIR")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
pytestmark = pytest.mark.skipif(
    not (APPLICATIONS and REFERENCE),
    reason="pinned MATLAB held-yaw output and BEMIO HDF5 not provided",
)


def _source_numerical_path(hydro, phases, *, end_time=250):
    """Replay source numerics from wave phases, never source motion or force.

    The published case has one moving yaw coordinate. WEC-Sim uses ode4 and
    extrapolates the last three major-step convolution loads with PCHIP at
    minor RK4 stages. This probe keeps that numerical route separate from
    the library's physical implicit-radiation default.
    """
    body = BodyClass(str(hydro))
    body.bodyNumber = 1
    body.bodyTotal = 2
    body.readH5file()
    data = body.hydroData
    components = pm_equal_energy_components(
        hydro, significant_height=2.5, peak_period=8,
        directions=[10], spreading=[1], phase=phases,
    )
    sampled = SampledPassiveYawExcitation.from_hydro_data(
        data, components, dt=0.005, end_time=end_time,
        ramp_time=100, rho=1000, g=9.81,
    )
    held = HeldPassiveYawExcitation(sampled, threshold=1)
    irf = data["hydro_coeffs"]["radiation_damping"]["impulse_response_fun"]
    kernel = CubicSpline(
        np.asarray(irf["t"]).ravel(),
        np.asarray(irf["K"])[5, 5] * 1000,
    )(np.arange(4001) * 0.01)
    added_mass = data["hydro_coeffs"]["added_mass"]["inf_freq"][5, 5] * 1000
    mass = 1.85e6 + added_mass
    dt = 0.01
    time = np.arange(round(end_time / dt) + 1) * dt
    yaw = np.zeros(len(time))
    speed = np.zeros(len(time))
    radiation = np.zeros(len(time))
    heading = np.zeros(len(time))
    radiation_times = np.array([-0.03, -0.02, -0.01])
    radiation_values = np.zeros(3)

    for step, at_time in enumerate(time):
        memory = min(step, 4000)
        history = speed[step - memory:step + 1][::-1]
        weights = np.ones(memory + 1)
        weights[0] = 0.5
        if step >= 4000:
            weights[-1] = 0.5
        # The source uses trapz over the entire 40 s IRF; early missing
        # velocity samples are zero, leaving the current sample half-weighted.
        radiation[step] = dt * np.dot(kernel[:memory + 1] * weights, history)
        radiation_times = np.r_[radiation_times[1:], at_time]
        radiation_values = np.r_[radiation_values[1:], radiation[step]]
        # The source's Memory block still exposes the preceding major-step
        # heading throughout this step's minor RK4 stages. Publish this
        # major-step selection only after all four stages have completed.
        current_heading = held._heading(yaw[step])
        heading[step] = current_heading
        if step == len(time) - 1:
            break

        extrapolated = PchipInterpolator(
            radiation_times, radiation_values, extrapolate=True,
        )

        def derivative(stage_time, angle, rate):
            torque = held(stage_time, np.array([angle]), np.array([rate]))[5]
            return rate, (torque - 120000 * rate
                          - float(extrapolated(stage_time))) / mass

        k1 = derivative(at_time, yaw[step], speed[step])
        k2 = derivative(at_time + dt / 2,
                        yaw[step] + dt * k1[0] / 2,
                        speed[step] + dt * k1[1] / 2)
        k3 = derivative(at_time + dt / 2,
                        yaw[step] + dt * k2[0] / 2,
                        speed[step] + dt * k2[1] / 2)
        k4 = derivative(at_time + dt,
                        yaw[step] + dt * k3[0],
                        speed[step] + dt * k3[1])
        yaw[step + 1] = yaw[step] + dt * (k1[0] + 2*k2[0] + 2*k3[0] + k4[0]) / 6
        speed[step + 1] = speed[step] + dt * (k1[1] + 2*k2[1] + 2*k3[1] + k4[1]) / 6
        held.last_heading = current_heading

    return time, yaw, speed, heading, sampled


@pytest.mark.parametrize("seed", (1, 2, 3))
def test_published_held_yaw_source_numerics_full_duration(seed):
    """Gate independent ode4 motion, heading events, and PTO work for 250 s."""
    reference = Path(REFERENCE)
    hydro = (Path(APPLICATIONS)
             / "_Common_Input_Files/OSWEC/hydroData/oswec.h5")
    prefix = "OSWEC_PASSIVE_YAW_IRR" + (f"_SEED{seed}" if seed > 1 else "")
    case = "PassiveYawRegression" + (f"Seed{seed}" if seed > 1 else "")
    components_csv = np.loadtxt(reference / f"{prefix}_components.csv", delimiter=",")
    source = np.loadtxt(
        reference / f"OSWEC_PASSIVE_YAW_IRR_{case}_body1.csv", delimiter=",",
    )
    pto = np.loadtxt(
        reference / f"OSWEC_PASSIVE_YAW_IRR_{case}_pto1.csv", delimiter=",",
    )
    wave = np.loadtxt(reference / f"{prefix}_wave.csv", delimiter=",")
    assert source.shape == pto.shape == (25001, 25)
    time, yaw, speed, heading, sampled = _source_numerical_path(
        hydro, components_csv[:, 3, None],
    )
    np.testing.assert_allclose(time, source[:, 0], rtol=0, atol=1e-10)
    np.testing.assert_allclose(sampled.elevation[::2], wave[:, 1],
                               rtol=0, atol=1e-11)
    np.testing.assert_allclose(pto[:, 17], -120000 * pto[:, 11],
                               rtol=0, atol=1e-5)
    # The source's major-step Memory latency is material: committing the
    # heading before the RK4 stages agrees early but separates by 250 s.
    yaw_gap = np.max(np.abs(yaw - source[:, 6]))
    speed_gap = np.max(np.abs(speed - source[:, 12]))
    torque_gap = np.max(np.abs(-120000 * speed - pto[:, 17]))
    assert yaw_gap < 1e-6
    assert speed_gap < 1e-7
    assert torque_gap < 0.01

    replay = HeldPassiveYawExcitation(sampled, threshold=1)
    source_heading = np.empty(len(time))
    for step, at_time in enumerate(time):
        replay.commit(at_time, np.array([source[step, 6]]))
        source_heading[step] = replay.last_heading
    source_events = np.flatnonzero(np.abs(np.diff(source_heading)) > 1e-5)
    python_events = np.flatnonzero(np.abs(np.diff(heading)) > 1e-5)
    np.testing.assert_array_equal(python_events, source_events)

    source_work = -np.trapezoid(pto[:, 23], time)
    python_work = 120000 * np.trapezoid(speed**2, time)
    work_gap = abs(python_work - source_work) / source_work
    assert source_work > 0 and np.isfinite((yaw, speed)).all()
    assert np.isfinite((yaw_gap, work_gap)).all()
    assert work_gap < 1e-6
    print(f"seed {seed}: {len(source_events)} heading updates paired; "
          f"250 s yaw {yaw_gap:.6g} rad, speed {speed_gap:.6g} rad/s, "
          f"PTO torque {torque_gap:.6g} N m, work {work_gap:.6g} relative")

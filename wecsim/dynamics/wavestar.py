"""Planar linkage used by the published WECCCOMP WaveStar model."""

from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np
from scipy.sparse import block_diag, csr_matrix

from wecsim.environment.irregularWave import IrregularComponents, synthesize_irregular_response
from wecsim.dynamics.rkdp87 import A as ODE8_A, B as ODE8_B, C as ODE8_C


@dataclass(frozen=True)
class WaveStarLinkage:
    """Rotate the float and arm about A and measure the B-to-C PTO stroke.

    Points are ``(x, z)`` world coordinates in the neutral pose. The PTO
    stroke is positive when the B-to-C distance shortens, as in the source.
    """

    pivot: tuple[float, float] = (-0.438, 0.302)
    pto_base: tuple[float, float] = (-0.438, 0.714)
    pto_arm: tuple[float, float] = (-0.6214398, 0.3816858)

    def point_position(self, neutral_point, angle):
        offset = np.asarray(neutral_point, dtype=float) - self.pivot
        cosine, sine = np.cos(angle), np.sin(angle)
        return np.asarray(self.pivot) + np.stack((
            cosine * offset[0] + sine * offset[1],
            -sine * offset[0] + cosine * offset[1],
        ), axis=-1)

    def point_velocity(self, neutral_point, angle, angular_speed):
        offset = np.asarray(neutral_point, dtype=float) - self.pivot
        cosine, sine = np.cos(angle), np.sin(angle)
        return np.asarray(angular_speed)[..., None] * np.stack((
            -sine * offset[0] + cosine * offset[1],
            -cosine * offset[0] - sine * offset[1],
        ), axis=-1)

    def pto_stroke(self, angle):
        initial_length = np.linalg.norm(
            np.asarray(self.pto_arm) - self.pto_base,
        )
        current_length = np.linalg.norm(
            self.point_position(self.pto_arm, angle) - self.pto_base,
            axis=-1,
        )
        return initial_length - current_length

    def pto_speed(self, angle, angular_speed):
        arm = self.point_position(self.pto_arm, angle) - self.pto_base
        arm_speed = self.point_velocity(self.pto_arm, angle, angular_speed)
        return -np.sum(arm * arm_speed, axis=-1) / np.linalg.norm(arm, axis=-1)


@dataclass(frozen=True)
class WaveStarResponse:
    time: np.ndarray
    wave_elevation: np.ndarray
    angle: np.ndarray
    angular_speed: np.ndarray
    float_position: np.ndarray
    float_velocity: np.ndarray
    arm_position: np.ndarray
    arm_velocity: np.ndarray
    pto_stroke: np.ndarray
    pto_speed: np.ndarray
    excitation_force: np.ndarray
    radiation_force: np.ndarray
    pto_force: np.ndarray


def run_wavestar_published(
    hydro_file: str | Path,
    components: IrregularComponents,
    *,
    dt: float = .01,
    end_time: float = 141.2,
    ramp_time: float = 7.06,
    rho: float = 1000.,
    g: float = 9.81,
    pto_controller=None,
    sensor_noise=None,
    sensor_dropout=None,
    fault_joint_friction=None,
    output_stride: int = 1,
    solver: str = "rk4",
) -> WaveStarResponse:
    """Integrate the published WECCCOMP WaveStar linkage.

    The pinned case uses a fitted state-space radiation model. Its active
    damping is positive, but the fit differs in magnitude from the BEM table.
    This runner uses that fit, with zero direct feedthrough as in the pinned
    WEC-Sim block, for comparison with the application. It does not change the
    general WEC radiation default. ``solver="wecsim_ode8"`` follows the
    pinned Simulink fixed-step solver and its delayed added-mass feedback;
    the default keeps the physical mass-matrix solve.
    """
    if (not np.isfinite([dt, end_time, ramp_time, rho, g]).all()
            or dt <= 0 or end_time <= 0 or ramp_time < 0
            or rho <= 0 or g <= 0):
        raise ValueError("WaveStar simulation settings are invalid")
    steps = round(end_time / dt)
    if not np.isclose(steps * dt, end_time, rtol=0, atol=1e-10):
        raise ValueError("WaveStar duration needs an integer number of steps")
    if (not isinstance(output_stride, int) or output_stride < 1
            or steps % output_stride):
        raise ValueError("WaveStar output_stride must divide the step count")
    if solver not in ("rk4", "wecsim_ode8"):
        raise ValueError("WaveStar solver must be rk4 or wecsim_ode8")
    if pto_controller is None:
        if sensor_noise is not None or sensor_dropout is not None:
            raise ValueError("WaveStar sensor inputs need a PTO controller")
    else:
        controller_step = getattr(pto_controller, "plant_dt", dt)
        if (not np.isfinite(controller_step)
                or not np.isclose(controller_step, dt, rtol=0, atol=1e-12)):
            raise ValueError("WaveStar PTO plant step must match the runner")
        if sensor_noise is None:
            sensor_noise = np.zeros(steps + 1)
        if sensor_dropout is None:
            sensor_dropout = np.zeros(steps + 1, dtype=bool)
        sensor_noise = np.asarray(sensor_noise, dtype=float)
        sensor_dropout = np.asarray(sensor_dropout, dtype=bool)
        if (sensor_noise.shape != (steps + 1,)
                or sensor_dropout.shape != (steps + 1,)
                or not np.isfinite(sensor_noise).all()):
            raise ValueError("WaveStar sensor records must cover every step")
    incident = synthesize_irregular_response(
        hydro_file, components, dt=dt / 2, end_time=end_time,
        ramp_time=ramp_time, rho=rho, g=g,
    )
    with h5py.File(hydro_file) as h5:
        center = np.asarray(h5["body1/properties/cg"]).ravel()
        buoyancy_center = np.asarray(h5["body1/properties/cb"]).ravel()
        volume = float(np.asarray(h5["body1/properties/disp_vol"]).item())
        if (center.shape != (3,) or buoyancy_center.shape != (3,)
                or not np.isfinite([*center, *buoyancy_center, volume]).all()
                or volume <= 0):
            raise ValueError("WaveStar hydrodynamic body properties are invalid")
        hydro = h5["body1/hydro_coeffs"]
        added_mass = np.asarray(hydro["added_mass/inf_freq"]) * rho
        restoring = (np.asarray(hydro["linear_restoring_stiffness"])
                     * rho * g)
        fit = hydro["radiation_damping/state_space"]
        raw_order = np.asarray(fit["it"])
        matrix_a = np.asarray(fit["A/all"])
        matrix_b = np.asarray(fit["B/all"])
        matrix_c = np.asarray(fit["C/all"]) * rho
    if (raw_order.shape != (6, 6) or not np.isfinite(raw_order).all()
            or not np.all(raw_order == np.floor(raw_order))
            or matrix_a.ndim != 4 or matrix_a.shape[:2] != (6, 6)
            or matrix_a.shape[2] != matrix_a.shape[3]
            or matrix_b.shape != (6, 6, matrix_a.shape[2], 1)
            or matrix_c.shape != (6, 6, 1, matrix_a.shape[2])
            or not all(np.isfinite(matrix).all()
                       for matrix in (matrix_a, matrix_b, matrix_c))
            or np.any(raw_order < 0)
            or np.any(raw_order > matrix_a.shape[2])):
        raise ValueError("WaveStar HDF5 needs a finite six-DOF radiation fit")
    order = raw_order.astype(int)
    if (added_mass.shape != (6, 6) or restoring.shape != (6, 6)
            or not np.isfinite(added_mass).all()
            or not np.isfinite(restoring).all()):
        raise ValueError("WaveStar HDF5 needs six-DOF hydrodynamics")

    blocks = []
    active = [(output, input_axis, int(order[output, input_axis]))
              for output in range(6) for input_axis in (0, 2, 4)
              if order[output, input_axis] > 0]
    count = sum(item[2] for item in active)
    input_matrix = np.zeros((count, 6))
    output_matrix = np.zeros((6, count))
    offset = 0
    for output, input_axis, size in active:
        blocks.append(csr_matrix(matrix_a[output, input_axis, :size, :size]))
        input_matrix[offset:offset + size, input_axis] = (
            matrix_b[output, input_axis, :size, 0]
        )
        output_matrix[output, offset:offset + size] = (
            matrix_c[output, input_axis, 0, :size]
        )
        offset += size
    state_matrix = block_diag(blocks, format="csr")

    linkage = WaveStarLinkage()
    auxiliary_base = np.asarray(linkage.pto_base)
    auxiliary_tip0 = np.asarray(linkage.pto_arm)
    auxiliary_rod0 = auxiliary_tip0 - auxiliary_base
    auxiliary_rod0_length_squared = auxiliary_rod0 @ auxiliary_rod0
    float_center = center[[0, 2]]
    arm_center = np.array([-.3301, .2551])
    float_mass, arm_mass = 3.075, 1.157
    float_inertia, arm_inertia = .001450, .0606
    rod_mass = motor_mass = rod_inertia = motor_inertia = .0001
    float_rigid = np.diag([float_mass] * 3 + [0., float_inertia, 0.])
    arm_rigid = np.diag([arm_mass] * 3 + [0., arm_inertia, 0.])
    source_ode8 = solver == "wecsim_ode8"
    if source_ode8:
        # WEC-Sim moves twice the translational A_inf trace into Simscape's
        # rigid mass. The remainder is fed back through a short delay.
        adjustment = 2 * np.trace(added_mass[:3, :3])
        float_effective = np.diag([float_mass + adjustment] * 3 +
                                  [0., float_inertia + added_mass[4, 4], 0.])
        delayed_added_mass = added_mass.copy()
        delayed_added_mass[:3, :3] -= adjustment * np.eye(3)
        for axis in range(3, 6):
            delayed_added_mass[axis, axis] = 0.
    else:
        float_effective = float_rigid + added_mass
    pitch_damping = np.diag([0., 0., 0., 0., 1.8, 0.])
    restoring_bias = np.zeros(6)
    restoring_bias[2] = (float_mass - rho * volume) * g
    restoring_bias[3] = rho * volume * g * (center[1] - buoyancy_center[1])
    restoring_bias[4] = rho * volume * g * (buoyancy_center[0] - center[0])
    arm_gravity = np.array([0., 0., -arm_mass * g, 0., 0., 0.])

    def body_motion(neutral_center, angle, speed):
        position = linkage.point_position(neutral_center, angle)
        tangent = linkage.point_velocity(neutral_center, angle, 1.)
        bias = -(position - linkage.pivot) * speed**2
        jacobian = np.array([tangent[0], 0., tangent[1], 0., 1., 0.])
        bias_acceleration = np.array([bias[0], 0., bias[1], 0., 0., 0.])
        return position, jacobian, bias_acceleration

    def auxiliary_dynamics(angle, speed):
        """Project the published BC rod and motor onto the arm coordinate.

        The rod translates to keep its C joint on the arm. Both light bodies
        rotate with the B-to-C line and contribute inertia and gravity.
        """
        tip = linkage.point_position(auxiliary_tip0, angle)
        tip_prime = linkage.point_velocity(auxiliary_tip0, angle, 1.)
        tip_second = -(tip - linkage.pivot)
        rod = tip - auxiliary_base
        length_squared = rod @ rod
        scale = np.sqrt(auxiliary_rod0_length_squared * length_squared)
        cosine = (auxiliary_rod0 @ rod) / scale
        sine = (auxiliary_rod0[1] * rod[0]
                - auxiliary_rod0[0] * rod[1]) / scale
        rotated_tip0 = np.array([
            cosine * auxiliary_tip0[0] + sine * auxiliary_tip0[1],
            -sine * auxiliary_tip0[0] + cosine * auxiliary_tip0[1],
        ])
        rotated_base = np.array([
            cosine * auxiliary_base[0] + sine * auxiliary_base[1],
            -sine * auxiliary_base[0] + cosine * auxiliary_base[1],
        ])
        tip_tangent = np.array([rotated_tip0[1], -rotated_tip0[0]])
        base_tangent = np.array([rotated_base[1], -rotated_base[0]])

        numerator = rod[1] * tip_prime[0] - rod[0] * tip_prime[1]
        theta_prime = numerator / length_squared
        numerator_prime = (rod[1] * tip_second[0]
                           - rod[0] * tip_second[1])
        theta_second = (numerator_prime * length_squared
                        - 2 * numerator * (rod @ tip_prime)) / length_squared**2
        rod_prime = tip_prime - theta_prime * tip_tangent
        motor_prime = -theta_prime * base_tangent
        rod_second = (tip_second - theta_second * tip_tangent
                      + theta_prime**2 * rotated_tip0)
        motor_second = (-theta_second * base_tangent
                        + theta_prime**2 * rotated_base)
        inertia = (rod_mass * (rod_prime @ rod_prime)
                   + motor_mass * (motor_prime @ motor_prime)
                   + (rod_inertia + motor_inertia) * theta_prime**2)
        inertia_prime = 2 * (
            rod_mass * (rod_prime @ rod_second)
            + motor_mass * (motor_prime @ motor_second)
            + (rod_inertia + motor_inertia) * theta_prime * theta_second
        )
        load = (.5 * inertia_prime * speed**2
                + g * (rod_mass * rod_prime[1]
                       + motor_mass * motor_prime[1]))
        return inertia, load

    def derivative(state, excitation, time, pto_force,
                   delayed_acceleration=None):
        angle, speed = state[:2]
        float_position, float_jacobian, float_bias = body_motion(
            float_center, angle, speed,
        )
        _, arm_jacobian, arm_bias = body_motion(arm_center, angle, speed)
        displacement = np.array([
            float_position[0] - float_center[0], 0.,
            float_position[1] - float_center[1], 0., angle, 0.,
        ])
        radiation = output_matrix @ state[2:]
        force = (excitation - radiation - restoring_bias
                 - restoring @ displacement
                 - pitch_damping @ (float_jacobian * speed))
        mass = (float_jacobian @ float_effective @ float_jacobian
                + arm_jacobian @ arm_rigid @ arm_jacobian)
        torque = (float_jacobian @ (force - float_effective @ float_bias)
                  + arm_jacobian @ (arm_gravity - arm_rigid @ arm_bias))
        auxiliary_mass, auxiliary_load = auxiliary_dynamics(angle, speed)
        mass += auxiliary_mass
        torque -= auxiliary_load
        if delayed_acceleration is not None:
            torque -= float_jacobian @ (delayed_added_mass @ delayed_acceleration)
        torque += pto_force * linkage.pto_speed(angle, 1.)
        if fault_joint_friction is not None:
            rod = (linkage.point_position(linkage.pto_arm, angle)
                   - linkage.pto_base)
            rod_tangent = linkage.point_velocity(linkage.pto_arm, angle, 1.)
            rod_angle_jacobian = ((rod[0] * rod_tangent[1]
                                   - rod[1] * rod_tangent[0]) / (rod @ rod))
            joint_jacobians = (-1 - rod_angle_jacobian,
                               -rod_angle_jacobian, 1.)
            for joint, joint_jacobian in enumerate(joint_jacobians):
                if joint == 2 and not (55 < time < 115):
                    continue
                torque += (fault_joint_friction.torque(joint_jacobian * speed)
                           * joint_jacobian)
        radiation_state = (state_matrix @ state[2:]
                           + input_matrix @ (float_jacobian * speed))
        return np.r_[speed, torque / mass, radiation_state]

    state = np.zeros(count + 2)
    sample_count = steps // output_stride + 1
    history = np.zeros((sample_count, count + 2))
    pto_forces = np.zeros(sample_count)
    body_acceleration_history = []
    # The NMPC motor holds torque; its axial force follows rod geometry.
    force_at_stroke = (None if pto_controller is None else
                       getattr(pto_controller, "force_at_stroke", None))
    for index in range(steps + 1):
        pto_force = (0. if pto_controller is None else pto_controller.step(
            linkage.pto_stroke(state[0]),
            noise=sensor_noise[index], dropout=sensor_dropout[index],
        ))
        if index % output_stride == 0:
            history[index // output_stride] = state
            pto_forces[index // output_stride] = pto_force
        if index == steps:
            break
        time = index * dt
        if source_ode8:
            stage = np.zeros((13, state.size))
            for number, fraction in enumerate(ODE8_C):
                trial = state + dt * (ODE8_A[number, :number] @ stage[:number])
                # The source delay uses the preceding accepted accelerations.
                # At the first stage of a major step, its logged value still
                # belongs to the preceding interval; the remaining stages
                # linearly extrapolate from the newest accepted acceleration.
                if fraction == 0:
                    delayed = (np.zeros(6) if index == 0 else
                               body_acceleration_history[0] if index == 1 else
                               2 * body_acceleration_history[index - 1]
                               - body_acceleration_history[index - 2])
                elif index == 0:
                    delayed = body_acceleration_history[0]
                else:
                    delayed = ((1 + fraction) * body_acceleration_history[index]
                               - fraction * body_acceleration_history[index - 1])
                sampled_index = 2 * index + 2 * fraction
                lower = int(np.floor(sampled_index))
                upper = min(lower + 1, len(incident.excitation_force) - 1)
                weight = sampled_index - lower
                excitation = ((1 - weight) * incident.excitation_force[lower]
                              + weight * incident.excitation_force[upper])
                stage_force = (pto_force if force_at_stroke is None else
                               force_at_stroke(linkage.pto_stroke(trial[0])))
                stage[number] = derivative(
                    trial, excitation, time + fraction * dt,
                    stage_force, delayed,
                )
                if number == 0:
                    _, jacobian, bias = body_motion(
                        float_center, state[0], state[1],
                    )
                    body_acceleration_history.append(
                        jacobian * stage[0, 1] + bias
                    )
            state += dt * (ODE8_B @ stage)
            continue
        first = derivative(state, incident.excitation_force[2 * index],
                           time, pto_force)
        second_state = state + dt * first / 2
        second_force = (pto_force if force_at_stroke is None else
                        force_at_stroke(linkage.pto_stroke(second_state[0])))
        second = derivative(second_state,
                            incident.excitation_force[2 * index + 1],
                            time + dt / 2, second_force)
        third_state = state + dt * second / 2
        third_force = (pto_force if force_at_stroke is None else
                       force_at_stroke(linkage.pto_stroke(third_state[0])))
        third = derivative(third_state,
                           incident.excitation_force[2 * index + 1],
                           time + dt / 2, third_force)
        fourth_state = state + dt * third
        fourth_force = (pto_force if force_at_stroke is None else
                        force_at_stroke(linkage.pto_stroke(fourth_state[0])))
        fourth = derivative(fourth_state,
                            incident.excitation_force[2 * index + 2],
                            time + dt, fourth_force)
        state += dt * (first + 2 * second + 2 * third + fourth) / 6

    angle = history[:, 0]
    speed = history[:, 1]

    def body_output(neutral_center):
        positions = np.zeros((sample_count, 6))
        velocities = np.zeros_like(positions)
        positions[:, [0, 2]] = linkage.point_position(neutral_center, angle)
        positions[:, 4] = angle
        velocities[:, [0, 2]] = linkage.point_velocity(
            neutral_center, angle, speed,
        )
        velocities[:, 4] = speed
        return positions, velocities

    float_position, float_velocity = body_output(float_center)
    arm_position, arm_velocity = body_output(arm_center)
    return WaveStarResponse(
        time=np.arange(sample_count) * (dt * output_stride),
        wave_elevation=incident.elevation[::2 * output_stride],
        angle=angle,
        angular_speed=speed,
        float_position=float_position,
        float_velocity=float_velocity,
        arm_position=arm_position,
        arm_velocity=arm_velocity,
        pto_stroke=linkage.pto_stroke(angle),
        pto_speed=linkage.pto_speed(angle, speed),
        excitation_force=incident.excitation_force[::2 * output_stride],
        radiation_force=history[:, 2:] @ output_matrix.T,
        pto_force=pto_forces,
    )

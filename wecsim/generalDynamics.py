"""Generalized-coordinate dynamics for supported linear-hydrodynamic devices.

The mechanical layout is supplied through body kinematics. Each body maps
generalized coordinates into its six WEC-Sim rigid-body coordinates, while
the engine assembles rigid inertia, cross-body added mass, restoring,
excitation, radiation, and PTO forces. Constant-frequency radiation uses RK4
or adaptive integration for continuous nonlinear forces; an impulse-response
kernel uses an implicit trapezoidal step.
"""

from dataclasses import dataclass
from typing import Callable

import numpy as np
from scipy.integrate import solve_ivp

from .directLinearGenerator import DirectLinearGenerator

@dataclass(frozen=True)
class BodyMotion:
    """Six-DOF displacement and its generalized-coordinate derivatives."""

    displacement: np.ndarray
    jacobian: np.ndarray
    bias_acceleration: np.ndarray  # J-dot(q, v) times v


@dataclass(frozen=True)
class DynamicBody:
    rigid_mass: np.ndarray
    added_mass: tuple[np.ndarray, ...]
    damping: tuple[np.ndarray, ...]
    restoring: np.ndarray
    static_force: np.ndarray
    reference_position: np.ndarray
    motion: Callable[[np.ndarray, np.ndarray], BodyMotion]
    excitation: Callable[[float], np.ndarray]
    radiation_kernel: np.ndarray | None = None  # (lag, 6, 6 * body_count)
    state_excitation: Callable[[float, np.ndarray, np.ndarray], np.ndarray] | None = None
    state_inertia: Callable[[float, np.ndarray, np.ndarray], tuple[np.ndarray, np.ndarray]] | None = None


@dataclass(frozen=True)
class DynamicsResponse:
    time: np.ndarray
    coordinate: np.ndarray
    speed: np.ndarray
    acceleration: np.ndarray
    body_position: np.ndarray
    body_velocity: np.ndarray
    controlled_pto_force: np.ndarray | None = None
    linear_generator_state: np.ndarray | None = None
    linear_generator_force: np.ndarray | None = None


class GeneralizedDynamics:
    """Integrate a constrained device described by independent coordinates.

    Body kinematics encode the constraints: the engine never integrates an
    unconstrained coordinate and then projects it back onto a joint. This
    avoids artificial constraint drift and accepts new layouts without
    changing the hydrodynamic force assembly.
    """

    def __init__(
        self,
        bodies: tuple[DynamicBody, ...],
        coordinate_count: int,
        *,
        pto_stiffness: np.ndarray | None = None,
        pto_damping: np.ndarray | None = None,
        pto_equilibrium: np.ndarray | None = None,
        pto_bias: np.ndarray | None = None,
        radiation_discretization: str = "trapezoid",
        added_mass_delay: float | None = None,
        nonlinear_force: Callable[[np.ndarray, np.ndarray], np.ndarray] | None = None,
        controlled_ptos: tuple = (),
        linear_generators: tuple = (),
    ):
        if not bodies or coordinate_count < 1:
            raise ValueError("a device needs bodies and independent coordinates")
        self.bodies = tuple(bodies)
        self.coordinate_count = coordinate_count
        n = coordinate_count
        self.pto_stiffness = (
            np.zeros((n, n)) if pto_stiffness is None
            else np.asarray(pto_stiffness, dtype=float)
        )
        self.pto_damping = (
            np.zeros((n, n)) if pto_damping is None
            else np.asarray(pto_damping, dtype=float)
        )
        self.pto_equilibrium = (
            np.zeros(n) if pto_equilibrium is None
            else np.asarray(pto_equilibrium, dtype=float)
        )
        self.pto_bias = (
            np.zeros(n) if pto_bias is None else np.asarray(pto_bias, dtype=float)
        )
        if radiation_discretization not in ("trapezoid", "fir"):
            raise ValueError("radiation_discretization must be trapezoid or fir")
        self.radiation_discretization = radiation_discretization
        if added_mass_delay is not None and (
            not np.isfinite(added_mass_delay) or added_mass_delay <= 0
            or radiation_discretization != "trapezoid"
            or any(body.radiation_kernel is None for body in bodies)
        ):
            raise ValueError("added-mass delay needs positive delay and convolution radiation")
        self.added_mass_delay = added_mass_delay
        if nonlinear_force is not None and not callable(nonlinear_force):
            raise TypeError("nonlinear_force must be callable")
        self.nonlinear_force = nonlinear_force
        self.controlled_ptos = tuple(controlled_ptos)
        self.linear_generators = tuple(linear_generators)
        for connection in self.controlled_ptos:
            if (connection.control is None
                    or np.shape(connection.stroke_jacobian) != (n,)
                    or not np.isfinite(connection.stroke_jacobian).all()):
                raise ValueError("controlled PTO needs a finite coordinate projection")
        for connection, generator in self.linear_generators:
            if (not isinstance(generator, DirectLinearGenerator)
                    or np.shape(connection.stroke_jacobian) != (n,)
                    or not np.isfinite(connection.stroke_jacobian).all()):
                raise ValueError("linear generator needs a finite PTO coordinate projection")
        if (self.pto_stiffness.shape != (n, n)
                or self.pto_damping.shape != (n, n)
                or self.pto_equilibrium.shape != (n,)
                or self.pto_bias.shape != (n,)
                or not np.isfinite(self.pto_stiffness).all()
                or not np.isfinite(self.pto_damping).all()
                or not np.isfinite(self.pto_equilibrium).all()
                or not np.isfinite(self.pto_bias).all()):
            raise ValueError("PTO matrices and force vectors must match the coordinate count")
        for body in self.bodies:
            if body.state_excitation is not None and not callable(body.state_excitation):
                raise TypeError("state_excitation must be callable")
            if body.state_inertia is not None and not callable(body.state_inertia):
                raise TypeError("state_inertia must be callable")
            if (len(body.added_mass) != len(bodies)
                    or len(body.damping) != len(bodies)):
                raise ValueError("each body needs one hydrodynamic block per body")
            for matrix in (body.rigid_mass, body.restoring,
                           *body.added_mass, *body.damping):
                if np.shape(matrix) != (6, 6) or not np.isfinite(matrix).all():
                    raise ValueError("body mass, restoring, and radiation blocks must be finite 6x6 matrices")
            for vector in (body.static_force, body.reference_position):
                if np.shape(vector) != (6,) or not np.isfinite(vector).all():
                    raise ValueError("body static force and reference position must be finite six-vectors")
            kernel = body.radiation_kernel
            if kernel is not None and (
                np.ndim(kernel) != 3 or np.shape(kernel)[0] < 1
                or np.shape(kernel)[1:] != (6, 6 * len(bodies))
                or not np.isfinite(kernel).all()
            ):
                raise ValueError("radiation kernel must have shape (lags, 6, 6 * bodies)")
        if (radiation_discretization == "fir"
                and any(body.radiation_kernel is None for body in self.bodies)):
            raise ValueError("FIR radiation needs a kernel for every body")
        if (self.controlled_ptos
                and any(body.radiation_kernel is not None for body in self.bodies)):
            raise ValueError("sampled PTO control currently needs constant radiation")
        if self.linear_generators and (
                self.controlled_ptos or self.added_mass_delay is not None
                or any(body.radiation_kernel is not None for body in self.bodies)):
            raise ValueError("linear generators currently need constant radiation and no sampled control")
        self.adjusted_rigid_mass = []
        self.applied_added_mass = []
        if added_mass_delay is not None:
            for index, body in enumerate(bodies):
                rigid = body.rigid_mass.copy()
                blocks = [block.copy() for block in body.added_mass]
                own = blocks[index]
                shift = np.zeros((6, 6))
                # WEC-Sim moves this part of A_inf into the Simscape body;
                # the remaining matrix acts on delayed body acceleration.
                shift[:3, :3] = 2 * np.trace(own[:3, :3]) * np.eye(3)
                for axis in range(3, 6):
                    shift[axis, axis] = own[axis, axis]
                    for other in range(axis + 1, 6):
                        shift[axis, other] = shift[other, axis] = own[axis, other]
                rigid += shift
                own -= shift
                self.adjusted_rigid_mass.append(rigid)
                self.applied_added_mass.append(tuple(blocks))
        self.delayed_body_acceleration = tuple(np.zeros(6) for _ in bodies)

    def acceleration(
        self,
        at_time: float,
        coordinate: np.ndarray,
        speed: np.ndarray,
        *,
        known_radiation: tuple[np.ndarray, ...] | None = None,
        dt: float | None = None,
        applied_force: np.ndarray | None = None,
    ) -> np.ndarray:
        """Assemble M(q) and generalized force, then solve M(q) q'' = F."""
        if any(body.radiation_kernel is not None for body in self.bodies):
            if known_radiation is None or dt is None or dt <= 0:
                raise ValueError("radiation-memory acceleration needs history and dt")
        if known_radiation is not None and len(known_radiation) != len(self.bodies):
            raise ValueError("radiation history needs one force vector per body")
        n = self.coordinate_count
        motions = tuple(body.motion(coordinate, speed) for body in self.bodies)
        mass = np.zeros((n, n))
        force = (-self.pto_stiffness @ (coordinate - self.pto_equilibrium)
                 - self.pto_damping @ speed + self.pto_bias)
        if self.nonlinear_force is not None:
            nonlinear = np.asarray(self.nonlinear_force(coordinate, speed),
                                   dtype=float)
            if nonlinear.shape != (n,) or not np.isfinite(nonlinear).all():
                raise ValueError("nonlinear force must be a finite generalized vector")
            force += nonlinear
        if applied_force is not None:
            added = np.asarray(applied_force, dtype=float)
            if added.shape != (n,) or not np.isfinite(added).all():
                raise ValueError("applied PTO force must be a finite generalized vector")
            force += added
        velocities = tuple(motion.jacobian @ speed for motion in motions)
        for i, body in enumerate(self.bodies):
            motion = motions[i]
            j = motion.jacobian
            rigid_mass = (body.rigid_mass if self.added_mass_delay is None
                          else self.adjusted_rigid_mass[i])
            excitation = (body.excitation(at_time)
                          if body.state_excitation is None else
                          np.asarray(body.state_excitation(
                              at_time, coordinate, speed), dtype=float))
            if excitation.shape != (6,) or not np.isfinite(excitation).all():
                raise ValueError("body excitation must be a finite six-vector")
            body_force = (body.static_force + excitation
                          - body.restoring @ motion.displacement
                          - rigid_mass @ motion.bias_acceleration)
            mass += j.T @ rigid_mass @ j
            if body.state_inertia is not None:
                extra_force, extra_mass = body.state_inertia(
                    at_time, coordinate, speed,
                )
                extra_force = np.asarray(extra_force, dtype=float)
                extra_mass = np.asarray(extra_mass, dtype=float)
                if (extra_force.shape != (6,) or extra_mass.shape != (6, 6)
                        or not np.isfinite(extra_force).all()
                        or not np.isfinite(extra_mass).all()):
                    raise ValueError("state inertia needs a finite force and 6x6 mass")
                body_force += extra_force - extra_mass @ motion.bias_acceleration
                mass += j.T @ extra_mass @ j
            for k, other in enumerate(motions):
                if self.added_mass_delay is None:
                    a = body.added_mass[k]
                    mass += j.T @ a @ other.jacobian
                    body_force -= a @ other.bias_acceleration
                else:
                    body_force -= (self.applied_added_mass[i][k]
                                   @ self.delayed_body_acceleration[k])
                body_force -= body.damping[k] @ velocities[k]
            if known_radiation is not None:
                body_force -= known_radiation[i]
                if (body.radiation_kernel is not None
                        and self.radiation_discretization == "trapezoid"):
                    current_velocity = np.concatenate(velocities)
                    body_force -= dt / 2 * (body.radiation_kernel[0] @ current_velocity)
            force += j.T @ body_force
        try:
            answer = np.linalg.solve(mass, force)
        except np.linalg.LinAlgError as exc:
            raise ValueError("the constrained device has a singular mass matrix") from exc
        if not np.isfinite(answer).all():
            raise ValueError("the dynamics produced nonfinite acceleration")
        return answer

    def integrate(
        self,
        *,
        dt: float,
        end_time: float,
        initial_coordinate: np.ndarray | None = None,
        initial_speed: np.ndarray | None = None,
        adaptive_regular: bool = False,
        applied_force_history: np.ndarray | None = None,
        step_force: Callable[[float, float, np.ndarray, np.ndarray], np.ndarray] | None = None,
    ) -> DynamicsResponse:
        """Integrate the body dynamics with optional sampled external loading.

        ``step_force`` receives the previous time, interval, and predicted
        next coordinate and speed. It advances once per radiation-memory
        step; its force acts during the ensuing implicit body solve. The
        initial external force is zero.
        """
        if not np.isfinite([dt, end_time]).all() or dt <= 0 or end_time < 0:
            raise ValueError("dt must be positive and end_time nonnegative")
        steps = round(end_time / dt)
        if not np.isclose(steps * dt, end_time, rtol=0, atol=1e-10):
            raise ValueError("end_time must be an integer multiple of dt")
        if self.added_mass_delay is not None and self.added_mass_delay >= dt:
            raise ValueError("added-mass delay must be shorter than the time step")
        if self.added_mass_delay is not None:
            self.delayed_body_acceleration = tuple(np.zeros(6) for _ in self.bodies)
        n = self.coordinate_count
        q = np.zeros((steps + 1, n))
        v = np.zeros_like(q)
        a = np.zeros_like(q)
        if initial_coordinate is not None:
            q[0] = np.asarray(initial_coordinate, dtype=float)
        if initial_speed is not None:
            v[0] = np.asarray(initial_speed, dtype=float)
        if not np.isfinite(q[0]).all() or not np.isfinite(v[0]).all():
            raise ValueError("initial state must be finite")
        time = np.arange(steps + 1) * dt
        memory = any(body.radiation_kernel is not None for body in self.bodies)
        if applied_force_history is not None:
            applied_force_history = np.asarray(applied_force_history, dtype=float)
            if (not memory or applied_force_history.shape != (steps + 1, n)
                    or not np.isfinite(applied_force_history).all()):
                raise ValueError("sampled applied force needs finite radiation-memory coordinate history")
        if step_force is not None:
            if not callable(step_force) or not memory or self.radiation_discretization != "trapezoid":
                raise ValueError("step_force needs callable and trapezoidal radiation memory")
            if applied_force_history is not None:
                raise ValueError("step_force and applied_force_history are mutually exclusive")
        if adaptive_regular and (memory or self.controlled_ptos
                                 or self.linear_generators
                                 or self.added_mass_delay is not None):
            raise ValueError("adaptive regular integration needs implicit mass, constant radiation, and no sampled control")
        controlled_force = None
        generator_state = None
        generator_force = None
        if adaptive_regular:
            self._integrate_adaptive_regular(time, q, v, a, dt)
        elif memory:
            if self.radiation_discretization == "fir":
                self._integrate_fir(time, q, v, a, dt)
            else:
                self._integrate_memory(time, q, v, a, dt, applied_force_history,
                                       step_force)
        elif self.controlled_ptos:
            controlled_force = self._integrate_controlled_regular(time, q, v, a, dt)
        elif self.linear_generators:
            generator_state, generator_force = self._integrate_linear_generator_regular(
                time, q, v, a, dt,
            )
        else:
            self._integrate_regular(time, q, v, a, dt)
        positions = np.zeros((steps + 1, len(self.bodies), 6))
        velocities = np.zeros_like(positions)
        for step in range(steps + 1):
            for i, body in enumerate(self.bodies):
                motion = body.motion(q[step], v[step])
                positions[step, i] = body.reference_position + motion.displacement
                velocities[step, i] = motion.jacobian @ v[step]
        return DynamicsResponse(time, q, v, a, positions, velocities,
                                controlled_force, generator_state,
                                generator_force)

    def _integrate_linear_generator_regular(self, time, q, v, a, dt):
        """Integrate body motion and continuous generator flux states together."""
        n = self.coordinate_count
        count = len(self.linear_generators)
        electrical = np.empty((len(time), count, 3))
        for index, (_, generator) in enumerate(self.linear_generators):
            electrical[0, index] = generator.initial_state()

        def forces_and_rates(speed, states):
            applied = np.zeros(n)
            rates = np.empty_like(states)
            forces = np.empty(count)
            for index, (connection, generator) in enumerate(self.linear_generators):
                stroke_speed = float(connection.stroke_jacobian @ speed)
                forces[index] = generator.force(stroke_speed, states[index])
                rates[index] = generator.state_rate(stroke_speed, states[index])
                applied += connection.stroke_jacobian * forces[index]
            return applied, rates, forces

        def derivative(at_time, state):
            coordinate = state[:n]
            speed = state[n:2*n]
            states = state[2*n:].reshape(count, 3)
            applied, rates, _ = forces_and_rates(speed, states)
            return np.r_[
                speed,
                self.acceleration(at_time, coordinate, speed,
                                  applied_force=applied),
                rates.ravel(),
            ]

        for step in range(len(time) - 1):
            at_time = time[step]
            state = np.r_[q[step], v[step], electrical[step].ravel()]
            k1 = derivative(at_time, state)
            k2 = derivative(at_time + dt / 2, state + dt * k1 / 2)
            k3 = derivative(at_time + dt / 2, state + dt * k2 / 2)
            k4 = derivative(at_time + dt, state + dt * k3)
            next_state = state + dt * (k1 + 2*k2 + 2*k3 + k4) / 6
            q[step + 1] = next_state[:n]
            v[step + 1] = next_state[n:2*n]
            electrical[step + 1] = next_state[2*n:].reshape(count, 3)
        force_history = np.empty((len(time), count))
        for step, at_time in enumerate(time):
            applied, _, forces = forces_and_rates(v[step], electrical[step])
            force_history[step] = forces
            a[step] = self.acceleration(
                at_time, q[step], v[step], applied_force=applied,
            )
        return electrical, force_history

    def _integrate_controlled_regular(self, time, q, v, a, dt):
        """Integrate regular-wave dynamics with sampled PTO controller memory."""
        n = self.coordinate_count
        states = [connection.control.initial_state()
                  for connection in self.controlled_ptos]
        force_history = np.zeros((len(time), len(states)))

        def controller_force(speed, memories):
            projected = np.zeros(n)
            forces = []
            next_memories = []
            for connection, memory in zip(self.controlled_ptos, memories):
                stroke_speed = float(connection.stroke_jacobian @ speed)
                force, next_memory = connection.control.sample(
                    stroke_speed, memory, dt,
                )
                projected += connection.stroke_jacobian * force
                forces.append(force)
                next_memories.append(next_memory)
            return projected, forces, next_memories

        for step, at_time in enumerate(time):
            applied, forces, next_states = controller_force(v[step], states)
            force_history[step] = forces
            a[step] = self.acceleration(
                at_time, q[step], v[step], applied_force=applied,
            )
            if step == len(time) - 1:
                break

            def derivative(stage_time, state):
                coordinate, speed = state[:n], state[n:]
                # Simulink's Memory blocks expose the previous major-step
                # controller state through intermediate RK4 evaluations.
                stage_force, _, _ = controller_force(speed, states)
                return np.concatenate((
                    speed,
                    self.acceleration(stage_time, coordinate, speed,
                                      applied_force=stage_force),
                ))

            state = np.concatenate((q[step], v[step]))
            k1 = derivative(at_time, state)
            k2 = derivative(at_time + dt / 2, state + dt * k1 / 2)
            k3 = derivative(at_time + dt / 2, state + dt * k2 / 2)
            k4 = derivative(at_time + dt, state + dt * k3)
            next_state = state + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6
            q[step + 1], v[step + 1] = next_state[:n], next_state[n:]
            states = next_states
        return force_history

    def _integrate_regular(self, time, q, v, a, dt):
        n = self.coordinate_count

        def derivative(at_time, state):
            coordinate, speed = state[:n], state[n:]
            return np.concatenate((
                speed, self.acceleration(at_time, coordinate, speed),
            ))

        for step in range(len(time) - 1):
            t = time[step]
            state = np.concatenate((q[step], v[step]))
            k1 = derivative(t, state)
            k2 = derivative(t + dt / 2, state + dt * k1 / 2)
            k3 = derivative(t + dt / 2, state + dt * k2 / 2)
            k4 = derivative(t + dt, state + dt * k3)
            next_state = state + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6
            q[step + 1], v[step + 1] = next_state[:n], next_state[n:]
        for step, t in enumerate(time):
            a[step] = self.acceleration(t, q[step], v[step])

    def _integrate_adaptive_regular(self, time, q, v, a, dt):
        """Resolve continuous nonlinear forces between requested output times."""
        n = self.coordinate_count

        def derivative(at_time, state):
            coordinate, speed = state[:n], state[n:]
            return np.concatenate((
                speed, self.acceleration(at_time, coordinate, speed),
            ))

        if len(time) > 1:
            solution = solve_ivp(
                derivative, (time[0], time[-1]),
                np.concatenate((q[0], v[0])), t_eval=time,
                method="RK45", rtol=1e-8, atol=1e-10, max_step=dt / 4,
            )
            if not solution.success or solution.y.shape != (2 * n, len(time)):
                raise RuntimeError("adaptive regular-wave integration failed")
            q[:] = solution.y[:n].T
            v[:] = solution.y[n:].T
        for step, at_time in enumerate(time):
            a[step] = self.acceleration(at_time, q[step], v[step])

    def _integrate_memory(self, time, q, v, a, dt, applied_force_history=None,
                          step_force=None):
        count = len(self.bodies)
        kernels = [body.radiation_kernel for body in self.bodies]
        if any(kernel is None for kernel in kernels):
            raise ValueError("all bodies must use the same radiation representation")
        velocity_history = np.zeros((len(time), count * 6))
        body_acceleration = (np.zeros((len(time), count, 6))
                             if self.added_mass_delay is not None else None)

        def save_body_state(step):
            for index, body in enumerate(self.bodies):
                motion = body.motion(q[step], v[step])
                velocity_history[step, 6 * index:6 * (index + 1)] = (
                    motion.jacobian @ v[step]
                )
                if body_acceleration is not None:
                    body_acceleration[step, index] = (
                        motion.jacobian @ a[step] + motion.bias_acceleration
                    )

        def known_radiation(step):
            result = []
            for kernel in kernels:
                memory = min(step, len(kernel) - 1)
                known = dt * np.einsum(
                    "tij,tj->i", kernel[1:memory + 1],
                    velocity_history[step - memory:step][::-1],
                )
                if len(kernel) > 1 and step >= len(kernel) - 1:
                    known -= dt / 2 * (
                        kernel[-1] @ velocity_history[step - memory]
                    )
                result.append(known)
            return tuple(result)

        def commit_sampled_excitation(step):
            for body in self.bodies:
                commit_state = getattr(body.state_excitation, "commit_state", None)
                if commit_state is not None:
                    commit_state(time[step], q[step], v[step])
                    continue
                commit = getattr(body.state_excitation, "commit", None)
                if commit is not None:
                    commit(time[step], q[step])

        zeros = tuple(np.zeros(6) for _ in self.bodies)
        a[0] = self.acceleration(time[0], q[0], v[0],
                                 known_radiation=zeros, dt=dt,
                                 applied_force=(None if applied_force_history is None
                                                else applied_force_history[0]))
        commit_sampled_excitation(0)
        save_body_state(0)
        for step in range(1, len(time)):
            if step_force is not None:
                # Native co-simulators may advance irreversibly. Predict the
                # next pose, then advance them exactly once before the body
                # solver's fixed-point iterations.
                predicted_speed = v[step - 1] + dt * a[step - 1]
                predicted_coordinate = (q[step - 1] + dt * v[step - 1]
                                        + dt * dt * a[step - 1] / 2)
                sampled_force = np.asarray(step_force(
                    time[step - 1], dt, predicted_coordinate,
                    predicted_speed), dtype=float)
                if (sampled_force.shape != (self.coordinate_count,)
                        or not np.isfinite(sampled_force).all()):
                    raise ValueError("step_force must return a finite generalized vector")
            else:
                sampled_force = (None if applied_force_history is None
                                 else applied_force_history[step])
            if body_acceleration is not None:
                if step == 1:
                    delayed = body_acceleration[0]
                else:
                    # Simulink's tiny Transport Delay extrapolates from the
                    # two preceding major-step acceleration samples.
                    factor = 1 - self.added_mass_delay / dt
                    delayed = (body_acceleration[step - 1] + factor
                               * (body_acceleration[step - 1]
                                  - body_acceleration[step - 2]))
                self.delayed_body_acceleration = tuple(delayed)
            known = known_radiation(step)
            trial_speed = v[step - 1].copy()
            # Mesh-pressure forces can make the trapezoidal fixed point
            # slower to converge near the instantaneous waterline.
            for _ in range(30):
                trial_coordinate = q[step - 1] + dt * (v[step - 1] + trial_speed) / 2
                trial_acceleration = self.acceleration(
                    time[step], trial_coordinate, trial_speed,
                    known_radiation=known, dt=dt,
                    applied_force=sampled_force,
                )
                next_speed = v[step - 1] + dt * (
                    a[step - 1] + trial_acceleration
                ) / 2
                if np.max(np.abs(next_speed - trial_speed)) < 1e-12:
                    trial_speed = next_speed
                    break
                trial_speed = next_speed
            else:
                raise RuntimeError(
                    f"radiation-memory step at t={time[step]:.6g} s did not "
                    "converge; check discontinuous forces or reduce dt"
                )
            v[step] = trial_speed
            q[step] = q[step - 1] + dt * (v[step - 1] + v[step]) / 2
            a[step] = self.acceleration(
                time[step], q[step], v[step], known_radiation=known, dt=dt,
                applied_force=sampled_force,
            )
            commit_sampled_excitation(step)
            save_body_state(step)

    def _integrate_fir(self, time, q, v, a, dt):
        """Sample the IRF at each step and hold its FIR force through RK4."""
        count = len(self.bodies)
        kernels = [body.radiation_kernel for body in self.bodies]
        velocity_history = np.zeros((len(time), count * 6))
        n = self.coordinate_count
        for step, at_time in enumerate(time):
            velocity_history[step] = np.concatenate([
                body.motion(q[step], v[step]).jacobian @ v[step]
                for body in self.bodies
            ])
            known = []
            for kernel in kernels:
                memory = min(step + 1, len(kernel))
                known.append(dt * np.einsum(
                    "tij,tj->i", kernel[:memory],
                    velocity_history[step - memory + 1:step + 1][::-1],
                ))
            known = tuple(known)
            a[step] = self.acceleration(
                at_time, q[step], v[step], known_radiation=known, dt=dt,
            )
            if step == len(time) - 1:
                break

            def derivative(stage_time, state):
                coordinate, speed = state[:n], state[n:]
                return np.concatenate((
                    speed,
                    self.acceleration(stage_time, coordinate, speed,
                                      known_radiation=known, dt=dt),
                ))

            state = np.concatenate((q[step], v[step]))
            k1 = derivative(at_time, state)
            k2 = derivative(at_time + dt / 2, state + dt * k1 / 2)
            k3 = derivative(at_time + dt / 2, state + dt * k2 / 2)
            k4 = derivative(at_time + dt, state + dt * k3)
            next_state = state + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6
            q[step + 1], v[step + 1] = next_state[:n], next_state[n:]

"""Three-phase direct linear generator from WEC-Sim's PTO-Sim block.

The source block uses a power-normalized ABC/dq transform, two stator flux
states, and an electrical angle. Its balanced resistive load can be evaluated
in dq coordinates; phase currents and voltages are reconstructed for output.
Force is signed as the force exerted by the generator on the PTO.
"""

from dataclasses import dataclass

import numpy as np

from wecsim.hydrodynamics.bodyClass import BodyClass


@dataclass(frozen=True)
class DirectLinearGeneratorSignals:
    force: float
    friction_force: float
    electromagnetic_force: float
    absorbed_power: float
    electrical_power: float
    current_d: float
    current_q: float
    phase_current: np.ndarray
    phase_voltage: np.ndarray


@dataclass(frozen=True)
class DirectLinearGenerator:
    """PTO-Sim's balanced, three-phase direct linear generator.

    ``load_resistance`` follows the source block's signed external-circuit
    convention: a passive load is negative. ``friction`` is likewise the
    signed force coefficient, so the published value is negative.
    """

    stator_resistance: float
    friction: float
    pole_pitch: float
    magnet_flux: float
    inductance: float
    load_resistance: float
    initial_angle: float = 0.0
    initial_flux_d: float | None = None
    initial_flux_q: float = 0.0

    def __post_init__(self):
        values = [self.stator_resistance, self.friction, self.pole_pitch,
                  self.magnet_flux, self.inductance, self.load_resistance,
                  self.initial_angle, self.initial_flux_q]
        if self.initial_flux_d is not None:
            values.append(self.initial_flux_d)
        if (not np.isfinite(values).all() or self.stator_resistance <= 0
                or self.friction > 0 or self.pole_pitch <= 0
                or self.magnet_flux <= 0 or self.inductance <= 0
                or self.load_resistance >= 0):
            raise ValueError("direct linear generator needs finite passive settings")

    def initial_state(self) -> np.ndarray:
        return np.array([
            self.magnet_flux if self.initial_flux_d is None
            else self.initial_flux_d,
            self.initial_flux_q,
            self.initial_angle,
        ], dtype=float)

    def _currents(self, state):
        state = np.asarray(state, dtype=float)
        if state.shape != (3,) or not np.isfinite(state).all():
            raise ValueError("generator state must be three finite values")
        return ((state[0] - self.magnet_flux) / self.inductance,
                state[1] / self.inductance)

    def state_rate(self, velocity: float, state: np.ndarray) -> np.ndarray:
        """Return d-flux, q-flux, and electrical-angle rates."""
        if not np.isfinite(velocity):
            raise ValueError("generator velocity must be finite")
        current_d, current_q = self._currents(state)
        omega = np.pi * velocity / self.pole_pitch
        resistance = self.load_resistance - self.stator_resistance
        return np.array([
            resistance * current_d + omega * state[1],
            resistance * current_q - omega * state[0],
            omega,
        ])

    def force(self, velocity: float, state: np.ndarray) -> float:
        """Force applied along the PTO stroke axis, without phase outputs."""
        if not np.isfinite(velocity):
            raise ValueError("generator velocity must be finite")
        _, current_q = self._currents(state)
        return self.friction * velocity + (
            np.pi / self.pole_pitch * self.magnet_flux * current_q
        )

    def signals(self, velocity: float, state: np.ndarray) -> DirectLinearGeneratorSignals:
        """Evaluate PTO force and balanced three-phase electrical output."""
        if not np.isfinite(velocity):
            raise ValueError("generator velocity must be finite")
        current_d, current_q = self._currents(state)
        angle = state[2] + np.array([0, -2 * np.pi / 3, 2 * np.pi / 3])
        phase_current = np.sqrt(2 / 3) * (
            current_d * np.cos(angle) - current_q * np.sin(angle)
        )
        phase_voltage = self.load_resistance * phase_current
        friction_force = self.friction * velocity
        electromagnetic_force = (
            np.pi / self.pole_pitch * self.magnet_flux * current_q
        )
        force = friction_force + electromagnetic_force
        return DirectLinearGeneratorSignals(
            force=force,
            friction_force=friction_force,
            electromagnetic_force=electromagnetic_force,
            absorbed_power=-force * velocity,
            electrical_power=-float(phase_voltage @ phase_current),
            current_d=current_d,
            current_q=current_q,
            phase_current=phase_current,
            phase_voltage=phase_voltage,
        )


@dataclass(frozen=True)
class RM3DirectGeneratorResponse:
    time: np.ndarray
    body_heave: np.ndarray  # time, float/spar, world z position
    body_heave_velocity: np.ndarray
    pto_stroke: np.ndarray
    pto_velocity: np.ndarray
    pto_force: np.ndarray
    absorbed_power: np.ndarray
    electrical_power: np.ndarray
    phase_current: np.ndarray  # time, a/b/c
    phase_voltage: np.ndarray
    friction_force: np.ndarray
    wave_elevation: np.ndarray
    flux_d: np.ndarray
    flux_q: np.ndarray
    electrical_angle: np.ndarray


def run_rm3_direct_linear_generator(
    hydro_file, *, generator: DirectLinearGenerator,
    height: float = 2.5, period: float = 8.0, ramp_time: float = 100.0,
    dt: float = 0.0005, output_dt: float = 0.01,
    end_time: float = 400.0, rho: float = 1000.0, g: float = 9.81,
) -> RM3DirectGeneratorResponse:
    """Run the published two-heave RM3 layout with a coupled generator.

    The source model has a world-to-spar translational constraint and a
    float-to-spar translational PTO, both along z. Other body motions are
    constrained. The BEM coefficients are read through production BodyClass.
    """
    settings = [height, period, ramp_time, dt, output_dt, end_time, rho, g]
    stride = round(output_dt / dt)
    steps = round(end_time / dt)
    if (not isinstance(generator, DirectLinearGenerator)
            or not np.isfinite(settings).all() or height < 0 or period <= 0
            or ramp_time < 0 or dt <= 0 or output_dt < dt or end_time <= 0
            or rho <= 0 or g <= 0 or stride < 1
            or not np.isclose(stride * dt, output_dt, rtol=0, atol=1e-10)
            or not np.isclose(steps * dt, end_time, rtol=0, atol=1e-9)
            or steps % stride):
        raise ValueError("RM3 direct-generator time and wave settings are invalid")
    omega = 2 * np.pi / period
    rigid_mass = np.empty(2)
    added_mass = np.empty(2)
    radiation = np.empty(2)
    restoring = np.empty(2)
    excitation_real = np.empty(2)
    excitation_imaginary = np.empty(2)
    equilibrium = np.empty(2)
    static_force = np.empty(2)
    for index in range(2):
        body = BodyClass(str(hydro_file))
        body.bodyNumber = index + 1
        body.bodyTotal = 2
        body.readH5file()
        if int(np.asarray(body.dof).item()) != 6:
            raise ValueError("the RM3 direct-generator model needs two six-DOF bodies")
        body.mass = "equilibrium"
        body.hydroStiffness = np.zeros((6, 6))
        body.viscDrag = {"Drag": np.zeros((6, 6)), "cd": np.zeros(6),
                         "characteristicArea": np.zeros(6)}
        body.linearDamping = np.zeros((6, 6))
        body.hydroForcePre(
            omega, [0], 1, np.array([0.0]), [], dt, rho, g, "regular",
            np.zeros((2, 1)), index + 1, 2, 0, 0, 0,
        )
        force = body.hydroForce
        rigid_mass[index] = float(np.asarray(body.mass).item())
        added_mass[index] = force["fAddedMass"][2, 2]
        radiation[index] = force["fDamping"][2, 2]
        restoring[index] = force["linearHydroRestCoef"][2, 2]
        excitation_real[index] = force["fExt"]["re"][2]
        excitation_imaginary[index] = force["fExt"]["im"][2]
        equilibrium[index] = np.asarray(body.cg).ravel()[2]
        static_force[index] = (
            rho * float(np.asarray(body.dispVol).item()) - rigid_mass[index]
        ) * g
    effective_mass = rigid_mass + added_mass
    if (not np.isfinite(np.r_[effective_mass, radiation, restoring,
                               excitation_real, excitation_imaginary,
                               equilibrium, static_force]).all()
            or np.any(effective_mass <= 0)):
        raise ValueError("RM3 heave hydrodynamics must be finite and nonsingular")

    def derivative(at_time, state):
        displacement = state[:2]
        velocity = state[2:4]
        electrical = state[4:]
        stroke_speed = velocity[0] - velocity[1]
        pto_force = generator.force(stroke_speed, electrical)
        ramp = (1.0 if ramp_time == 0 or at_time >= ramp_time else
                (1 - np.cos(np.pi * at_time / ramp_time)) / 2)
        excitation = height / 2 * ramp * (
            excitation_real * np.cos(omega * at_time)
            - excitation_imaginary * np.sin(omega * at_time)
        )
        acceleration = (
            excitation + static_force - restoring * displacement
            - radiation * velocity + pto_force * np.array([1.0, -1.0])
        ) / effective_mass
        return np.r_[velocity, acceleration,
                     generator.state_rate(stroke_speed, electrical)]

    time = np.arange(steps // stride + 1) * output_dt
    saved = np.empty((len(time), 7))
    state = np.r_[np.zeros(4), generator.initial_state()]
    saved[0] = state
    for step in range(steps):
        at_time = step * dt
        k1 = derivative(at_time, state)
        k2 = derivative(at_time + dt / 2, state + dt * k1 / 2)
        k3 = derivative(at_time + dt / 2, state + dt * k2 / 2)
        k4 = derivative(at_time + dt, state + dt * k3)
        state += dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6
        if (step + 1) % stride == 0:
            saved[(step + 1) // stride] = state
    displacement = saved[:, :2]
    velocity = saved[:, 2:4]
    stroke = displacement[:, 0] - displacement[:, 1]
    stroke_speed = velocity[:, 0] - velocity[:, 1]
    signals = [generator.signals(speed, electrical)
               for speed, electrical in zip(stroke_speed, saved[:, 4:])]
    ramp = np.ones_like(time)
    if ramp_time:
        early = time < ramp_time
        ramp[early] = (1 - np.cos(np.pi * time[early] / ramp_time)) / 2
    return RM3DirectGeneratorResponse(
        time=time,
        body_heave=displacement + equilibrium,
        body_heave_velocity=velocity,
        pto_stroke=stroke,
        pto_velocity=stroke_speed,
        pto_force=np.array([item.force for item in signals]),
        absorbed_power=np.array([item.absorbed_power for item in signals]),
        electrical_power=np.array([item.electrical_power for item in signals]),
        phase_current=np.array([item.phase_current for item in signals]),
        phase_voltage=np.array([item.phase_voltage for item in signals]),
        friction_force=np.array([item.friction_force for item in signals]),
        wave_elevation=height / 2 * ramp * np.cos(omega * time),
        flux_d=saved[:, 4], flux_q=saved[:, 5], electrical_angle=saved[:, 6],
    )

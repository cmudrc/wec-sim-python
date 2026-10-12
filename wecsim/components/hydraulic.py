"""Compressible hydraulic cylinder equations used by WEC-Sim PTO-Sim.

The components below can be assembled into a rectified PTO network for a
coupled WEC trajectory.
"""

from dataclasses import dataclass

import numpy as np

from wecsim.components.electricGenerator import DiscretePILoadController, EquivalentCircuitGenerator


@dataclass(frozen=True)
class CompressibleCylinder:
    area_a: float
    area_b: float
    bulk_modulus: float
    stroke: float
    offset: float

    def __post_init__(self):
        values = np.asarray((self.area_a, self.area_b, self.bulk_modulus,
                             self.stroke, self.offset), dtype=float)
        if not np.isfinite(values).all() or np.any(values <= 0):
            raise ValueError("cylinder areas, bulk modulus, stroke, and offset must be positive")

    def force(self, pressure_a, pressure_b):
        """Force in the pinned PTO-Sim output sign convention, in newtons."""
        return (self.area_b * np.asarray(pressure_b, dtype=float)
                - self.area_a * np.asarray(pressure_a, dtype=float))

    def pressure_rates(self, position, velocity, flow_a, flow_b):
        """Return A/B pressure rates for position, speed, and port flows.

        The volume floors match the pinned Simulink cylinder blocks. Inputs
        are SI metres, m/s, and m³/s; outputs are Pa/s.
        """
        position = np.asarray(position, dtype=float)
        velocity = np.asarray(velocity, dtype=float)
        flow_a = np.asarray(flow_a, dtype=float)
        flow_b = np.asarray(flow_b, dtype=float)
        length_a = np.maximum(self.offset - position, 1e-10)
        volume_b = np.maximum((position - self.offset + self.stroke)
                              * self.area_b, 1e-10)
        rate_a = (self.bulk_modulus * (self.area_a * velocity + flow_a)
                  / (self.area_a * length_a))
        rate_b = (self.bulk_modulus * (flow_b - self.area_b * velocity)
                  / volume_b)
        return rate_a, rate_b


@dataclass(frozen=True)
class IdealDoubleActingCylinder:
    """Incompressible two-chamber cylinder away from its travel stops.

    Rod speed is positive when chamber A expands. Port B flow is reported in
    the source block's outflow convention. Network pressures are inputs; this
    component does not solve the connected hydraulic circuit.
    """

    area_a: float
    area_b: float

    def __post_init__(self):
        areas = np.asarray((self.area_a, self.area_b), dtype=float)
        if not np.isfinite(areas).all() or np.any(areas <= 0):
            raise ValueError("cylinder piston areas must be finite and positive")

    def force(self, pressure_a, pressure_b):
        """Force exerted by the cylinder rod, in newtons."""
        return (self.area_a * np.asarray(pressure_a, dtype=float)
                - self.area_b * np.asarray(pressure_b, dtype=float))

    def port_flows(self, rod_speed):
        """A-port inflow and B-port outflow, in m³/s."""
        speed = np.asarray(rod_speed, dtype=float)
        return self.area_a * speed, self.area_b * speed


@dataclass(frozen=True)
class RectifyingCheckValve:
    """Four-port PTO-Sim valve with smooth pressure-dependent check openings."""

    discharge_coefficient: float
    area_max: float
    area_min: float
    pressure_max: float
    pressure_min: float
    density: float
    switch_gain: float
    opening_gain: float

    def __post_init__(self):
        values = np.asarray((self.discharge_coefficient, self.area_max,
                             self.area_min, self.pressure_max,
                             self.pressure_min, self.density,
                             self.switch_gain, self.opening_gain), dtype=float)
        if not np.isfinite(values).all() or (values[[0, 1, 5, 6, 7]] <= 0).any():
            raise ValueError("valve coefficient, maximum area, density, and gains must be positive")
        if (self.area_min < 0 or self.area_min > self.area_max
                or self.pressure_max <= self.pressure_min):
            raise ValueError("valve area and pressure bounds are invalid")

    def check_flow(self, pressure_drop):
        """Nonnegative flow through one source check valve, in m³/s."""
        drop = np.asarray(pressure_drop, dtype=float)
        midpoint = (self.pressure_max + self.pressure_min) / 2
        area = self.area_min + (self.area_max - self.area_min) / 2 * (
            1 + np.tanh(self.opening_gain * (drop - midpoint)))
        speed = np.sqrt(2 * drop * np.tanh(self.switch_gain * drop)
                        / self.density)
        return self.discharge_coefficient * area * speed

    def flows(self, pressure_a, pressure_b, pressure_high, pressure_low):
        """Return signed flows at cylinder A/B and high/low pressure ports."""
        a, b, high, low = np.broadcast_arrays(
            np.asarray(pressure_a, dtype=float),
            np.asarray(pressure_b, dtype=float),
            np.asarray(pressure_high, dtype=float),
            np.asarray(pressure_low, dtype=float),
        )
        a_to_high = self.check_flow(a - high)
        b_to_high = self.check_flow(b - high)
        low_to_b = self.check_flow(low - b)
        low_to_a = self.check_flow(low - a)
        return (low_to_a - a_to_high, low_to_b - b_to_high,
                a_to_high + b_to_high, -low_to_b - low_to_a)


@dataclass(frozen=True)
class GasChargedAccumulator:
    """Polytropic gas pressure versus liquid volume and inlet flow.

    Optional atmospheric pressure and finite hard stops describe the legacy
    Simscape Fluids accumulator. The defaults retain the PTO-Sim accumulator
    law, which uses gauge pressure and no volume hard stops.
    """

    initial_gas_volume: float
    precharge_pressure: float
    exponent: float = 1.4
    atmospheric_pressure: float = 0.0
    dead_gas_volume: float = 0.0
    hard_stop_stiffness: float = 0.0
    hard_stop_damping: float = 0.0

    def __post_init__(self):
        values = np.asarray((self.initial_gas_volume, self.precharge_pressure,
                             self.exponent, self.atmospheric_pressure,
                             self.dead_gas_volume, self.hard_stop_stiffness,
                             self.hard_stop_damping), dtype=float)
        if (not np.isfinite(values).all() or (values[:3] <= 0).any()
                or (values[3:] < 0).any()
                or self.dead_gas_volume >= self.initial_gas_volume):
            raise ValueError(
                "accumulator volume, precharge, and exponent must be positive; "
                "atmospheric pressure and hard stops must be nonnegative"
            )

    def pressure(self, liquid_inflow_volume, flow_rate=0.0):
        """Gauge pressure in Pa; positive liquid volume compresses gas."""
        liquid_volume, flow = np.broadcast_arrays(
            np.asarray(liquid_inflow_volume, dtype=float),
            np.asarray(flow_rate, dtype=float),
        )
        if not np.isfinite(liquid_volume).all() or not np.isfinite(flow).all():
            raise ValueError("accumulator volume and flow must be finite")
        gas_volume = self.initial_gas_volume - liquid_volume
        if np.any(gas_volume <= 0):
            raise ValueError("accumulator gas volume must remain positive")
        gas_pressure = ((self.precharge_pressure + self.atmospheric_pressure)
                        * (self.initial_gas_volume / gas_volume) ** self.exponent
                        - self.atmospheric_pressure)
        capacity = self.initial_gas_volume - self.dead_gas_volume
        below = np.minimum(liquid_volume, 0)
        above = np.maximum(liquid_volume - capacity, 0)
        contact_pressure = (
            (self.hard_stop_stiffness
             - self.hard_stop_damping * np.minimum(flow, 0)) * below
            + (self.hard_stop_stiffness
               + self.hard_stop_damping * np.maximum(flow, 0)) * above
        )
        return gas_pressure + contact_pressure


@dataclass(frozen=True)
class ConstantEfficiencyHydraulicMotor:
    """PTO-Sim hydraulic motor with constant volumetric/mechanical efficiency."""

    displacement_m3_per_rev: float
    volumetric_efficiency: float
    mechanical_efficiency: float

    def __post_init__(self):
        values = np.asarray((self.displacement_m3_per_rev,
                             self.volumetric_efficiency,
                             self.mechanical_efficiency), dtype=float)
        if not np.isfinite(values).all() or (values <= 0).any():
            raise ValueError("motor displacement and efficiencies must be positive")
        if self.volumetric_efficiency > 1 or self.mechanical_efficiency > 1:
            raise ValueError("motor efficiencies cannot exceed one")

    def torque(self, pressure_drop):
        """Shaft torque in N m from high-minus-low pressure, in Pa."""
        return (self.mechanical_efficiency * self.displacement_m3_per_rev
                * np.asarray(pressure_drop, dtype=float) / (2 * np.pi))

    def flow(self, angular_speed):
        """Hydraulic inlet flow in m³/s for shaft speed in rad/s."""
        return (self.displacement_m3_per_rev
                * np.asarray(angular_speed, dtype=float)
                / (2 * np.pi * self.volumetric_efficiency))


@dataclass(frozen=True)
class RectifiedHydraulicPTO:
    """RM3-style cylinder, valve, accumulator, motor, and generator network.

    State order is chamber pressures A/B (Pa), accumulated inlet volumes at
    high/low pressure (m³), shaft speed (rad/s), current (A), and PI integral
    resistance (ohm). The source uses forward Euler at each simulation step.
    """

    cylinder: CompressibleCylinder
    valve: RectifyingCheckValve
    high_accumulator: GasChargedAccumulator
    low_accumulator: GasChargedAccumulator
    motor: ConstantEfficiencyHydraulicMotor
    generator: EquivalentCircuitGenerator
    load_controller: DiscretePILoadController
    initial_pressure_a: float
    initial_pressure_b: float
    initial_shaft_speed: float = 0.0
    initial_current: float = 0.0
    initial_load_integral: float = 0.0

    def __post_init__(self):
        expected = (
            (self.cylinder, CompressibleCylinder),
            (self.valve, RectifyingCheckValve),
            (self.high_accumulator, GasChargedAccumulator),
            (self.low_accumulator, GasChargedAccumulator),
            (self.motor, ConstantEfficiencyHydraulicMotor),
            (self.generator, EquivalentCircuitGenerator),
            (self.load_controller, DiscretePILoadController),
        )
        if any(not isinstance(value, cls) for value, cls in expected):
            raise TypeError("hydraulic PTO components have invalid types")
        values = np.asarray((self.initial_pressure_a, self.initial_pressure_b,
                             self.initial_shaft_speed, self.initial_current,
                             self.initial_load_integral), dtype=float)
        if (not np.isfinite(values).all() or (values[:2] <= 0).any()):
            raise ValueError("hydraulic PTO initial state must be finite with positive pressures")

    def initial_state(self):
        return np.array([self.initial_pressure_a, self.initial_pressure_b,
                         0., 0., self.initial_shaft_speed,
                         self.initial_current, self.initial_load_integral])

    def force(self, state):
        values = np.asarray(state, dtype=float)
        if values.ndim == 0 or values.shape[-1] != 7:
            raise ValueError("hydraulic PTO state must have seven entries")
        return self.cylinder.force(values[..., 0], values[..., 1])

    def state_rate(self, stroke, stroke_speed, state):
        values = np.asarray(state, dtype=float)
        if values.shape != (7,) or not np.isfinite(values).all():
            raise ValueError("hydraulic PTO state must have seven finite entries")
        p_a, p_b, volume_high, volume_low, omega, current, integral = values
        p_high = self.high_accumulator.pressure(volume_high)
        p_low = self.low_accumulator.pressure(volume_low)
        flow_a, flow_b, flow_high, flow_low = self.valve.flows(
            p_a, p_b, p_high, p_low)
        motor_flow = self.motor.flow(omega)
        motor_torque = self.motor.torque(p_high - p_low)
        shaft_rpm = omega * 60 / (2 * np.pi)
        voltage = self.load_controller.voltage(shaft_rpm, current, integral)
        rate_a, rate_b = self.cylinder.pressure_rates(
            stroke, stroke_speed, flow_a, flow_b)
        return np.array([
            rate_a, rate_b, flow_high - motor_flow, flow_low + motor_flow,
            self.generator.speed_rate(omega, current, motor_torque),
            self.generator.current_rate(omega, current, voltage),
            self.load_controller.integral_rate(shaft_rpm),
        ], dtype=float)

"""Hydraulic components in the published OSWEC reverse-osmosis application."""

from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq

from wecsim.components.hydraulic import GasChargedAccumulator, IdealDoubleActingCylinder


@dataclass(frozen=True)
class ReverseOsmosisMembrane:
    """Osmotic relief valve and hydraulic resistance in series.

    Pressures are gauge pressure relative to the permeate outlet. The valve
    follows the published hydraulic relief-block opening and laminar
    transition equations. ``permeate_flow`` solves the pressure split across
    the valve and linear resistance; it does not determine network pressure.
    """

    resistance: float
    osmotic_pressure: float
    regulation_range: float
    max_area: float
    leakage_area: float
    discharge_coefficient: float
    fluid_density: float
    laminar_pressure_ratio: float = .999
    atmospheric_pressure: float = 101_325.

    def __post_init__(self):
        values = np.asarray((
            self.resistance, self.osmotic_pressure, self.regulation_range,
            self.max_area, self.leakage_area, self.discharge_coefficient,
            self.fluid_density, self.laminar_pressure_ratio,
            self.atmospheric_pressure,
        ), dtype=float)
        if (not np.isfinite(values).all() or np.any(values[:7] <= 0)
                or self.leakage_area > self.max_area
                or not 0 < self.laminar_pressure_ratio < 1
                or self.atmospheric_pressure <= 0):
            raise ValueError("reverse-osmosis valve and resistance settings are invalid")

    def valve_flow(self, pressure_drop: float) -> float:
        """Signed flow through the pressure-controlled osmotic valve."""
        dp = float(pressure_drop)
        if not np.isfinite(dp):
            raise ValueError("valve pressure drop must be finite")
        opening = np.clip(
            (dp - self.osmotic_pressure) / self.regulation_range, 0, 1,
        )
        area = self.leakage_area + (self.max_area - self.leakage_area) * opening
        critical = ((self.atmospheric_pressure + dp / 2)
                    * (1 - self.laminar_pressure_ratio))
        return (self.discharge_coefficient * area
                * np.sqrt(2 / self.fluid_density) * dp
                / (dp * dp + critical * critical) ** .25)

    def permeate_flow(self, inlet_pressure):
        """Solve permeate flow from the measured or modeled inlet pressure."""
        pressure = np.asarray(inlet_pressure, dtype=float)
        if not np.isfinite(pressure).all() or np.any(pressure < 0):
            raise ValueError("membrane inlet pressure must be finite and nonnegative")
        flows = np.empty_like(pressure)
        for index, value in np.ndenumerate(pressure):
            if value == 0:
                flows[index] = 0
                continue
            drop = brentq(
                lambda dp: dp + self.resistance * self.valve_flow(dp) - value,
                0, float(value), xtol=1e-8,
            )
            flows[index] = self.valve_flow(drop)
        return float(flows) if flows.ndim == 0 else flows


@dataclass(frozen=True)
class DynamicPressureReliefValve:
    """Pressure relief orifice with the published first-order opening lag."""

    set_pressure: float
    regulation_range: float
    max_area: float
    leakage_area: float
    discharge_coefficient: float
    fluid_density: float
    opening_time_constant: float
    laminar_pressure_ratio: float = .999
    atmospheric_pressure: float = 101_325.

    def __post_init__(self):
        values = np.asarray((
            self.set_pressure, self.regulation_range, self.max_area,
            self.leakage_area, self.discharge_coefficient, self.fluid_density,
            self.opening_time_constant, self.laminar_pressure_ratio,
            self.atmospheric_pressure,
        ), dtype=float)
        if (not np.isfinite(values).all() or np.any(values <= 0)
                or self.leakage_area > self.max_area
                or not 0 < self.laminar_pressure_ratio < 1):
            raise ValueError("pressure relief settings are invalid")

    def next_area(self, pressure: float, previous_area: float, dt: float) -> float:
        """Backward Euler update of the pressure-controlled orifice area."""
        opening = np.clip((pressure - self.set_pressure)
                          / self.regulation_range, 0, 1)
        target = self.leakage_area + (self.max_area - self.leakage_area) * opening
        ratio = dt / self.opening_time_constant
        return float((previous_area + ratio * target) / (1 + ratio))

    def flow(self, pressure_drop: float, area: float) -> float:
        """Signed source orifice flow at the current dynamic opening."""
        dp = float(pressure_drop)
        critical = ((self.atmospheric_pressure + dp / 2)
                    * (1 - self.laminar_pressure_ratio))
        return float(self.discharge_coefficient * area
                     * np.sqrt(2 / self.fluid_density) * dp
                     / (dp * dp + critical * critical) ** .25)


@dataclass(frozen=True)
class ReverseOsmosisHydraulicState:
    """Pressure, liquid volume, valve opening, and flows after one network step."""

    pressure: float
    liquid_volume: float
    relief_area: float
    cylinder_feed_flow: float
    permeate_flow: float
    brine_flow: float
    recovered_feed_flow: float
    accumulator_flow: float
    relief_flow: float


@dataclass(frozen=True)
class ReverseOsmosisHydraulicNetwork:
    """Published OSWEC desalination high-pressure hydraulic junction.

    The ideal cylinder's two check-valve paths supply ``area * abs(speed)``.
    A motor driving a smaller pump recovers part of the brine flow. This
    network predicts high-side pressure and flows from prescribed rod speed;
    it does not determine cylinder chamber pressure or WEC motion.
    """

    accumulator: GasChargedAccumulator
    membrane: ReverseOsmosisMembrane
    relief: DynamicPressureReliefValve
    cylinder_area: float
    pump_to_motor_displacement_ratio: float
    pump_outlet_resistance: float

    def __post_init__(self):
        values = np.asarray((self.cylinder_area,
                             self.pump_to_motor_displacement_ratio,
                             self.pump_outlet_resistance), dtype=float)
        if (not isinstance(self.accumulator, GasChargedAccumulator)
                or not isinstance(self.membrane, ReverseOsmosisMembrane)
                or not isinstance(self.relief, DynamicPressureReliefValve)
                or not np.isfinite(values).all() or np.any(values <= 0)
                or self.pump_to_motor_displacement_ratio >= 1):
            raise ValueError("reverse-osmosis network settings are invalid")

    def initial_state(self, pressure: float = 0.) -> ReverseOsmosisHydraulicState:
        """Find accumulator liquid volume for an initial gauge pressure."""
        if not np.isfinite(pressure) or pressure < 0:
            raise ValueError("initial pressure must be finite and nonnegative")
        lower = -self.accumulator.initial_gas_volume
        for _ in range(32):
            if self.accumulator.pressure(lower, 0) <= pressure:
                break
            lower *= 2
        else:
            raise ValueError("could not bracket initial accumulator volume")
        upper = np.nextafter(self.accumulator.initial_gas_volume, 0.)
        volume = brentq(
            lambda value: self.accumulator.pressure(value, 0) - pressure,
            lower, upper,
        )
        return ReverseOsmosisHydraulicState(
            float(pressure), volume, self.relief.leakage_area,
            0., 0., 0., 0., 0., 0.,
        )

    def step(self, rod_speed: float, previous: ReverseOsmosisHydraulicState,
             dt: float) -> ReverseOsmosisHydraulicState:
        """Advance high-side pressure and flows from prescribed rod speed."""
        if not np.isfinite(rod_speed):
            raise ValueError("rod speed must be finite")
        return self.step_from_feed_flow(self.cylinder_area * abs(rod_speed),
                                        previous, dt)

    def step_from_feed_flow(self, feed_flow: float,
                            previous: ReverseOsmosisHydraulicState,
                            dt: float) -> ReverseOsmosisHydraulicState:
        """Advance from supplied feed flow, useful with measured cylinder flow."""
        if not np.isfinite(feed_flow) or not np.isfinite(dt) or dt <= 0:
            raise ValueError("feed flow must be finite and time step positive")
        ratio = self.pump_to_motor_displacement_ratio

        def flows(pressure):
            brine = (pressure * (1 - ratio)
                     / (ratio * ratio * self.pump_outlet_resistance))
            recovered = ratio * brine
            permeate = self.membrane.permeate_flow(pressure)
            area = self.relief.next_area(pressure, previous.relief_area, dt)
            relief = self.relief.flow(pressure, area)
            accumulator = feed_flow + recovered - brine - permeate - relief
            return accumulator, area, permeate, brine, recovered, relief

        def residual(pressure):
            accumulator, *_ = flows(pressure)
            return (self.accumulator.pressure(
                previous.liquid_volume + dt * accumulator, accumulator,
            ) - pressure)

        if residual(0) < 0:
            raise ValueError("network pressure would be negative")
        upper = max(1e6, 2 * previous.pressure)
        for _ in range(32):
            if residual(upper) <= 0:
                break
            upper *= 2
        else:
            raise ValueError("could not bracket nonnegative network pressure")
        pressure = brentq(residual, 0, upper, xtol=1e-6)
        accumulator, area, permeate, brine, recovered, relief = flows(pressure)
        return ReverseOsmosisHydraulicState(
            pressure, previous.liquid_volume + dt * accumulator, area,
            feed_flow, permeate, brine, recovered, accumulator, relief,
        )


@dataclass(frozen=True)
class FourValveRectifiedCylinder:
    """Two chamber pressures from prescribed rod speed and high/low lines.

    Each chamber connects to both lines through a two-way directional valve.
    Positive rod speed opens B-to-high and A-to-low; negative speed reverses
    the paths. This solves the documented passive orifice law, without the
    alternating pressure artifact in the pinned legacy Simscape trace.
    """

    cylinder: IdealDoubleActingCylinder
    max_area: float
    leakage_area: float
    discharge_coefficient: float
    fluid_density: float
    laminar_pressure_ratio: float = .999
    atmospheric_pressure: float = 101_325.

    def __post_init__(self):
        values = np.asarray((
            self.max_area, self.leakage_area, self.discharge_coefficient,
            self.fluid_density, self.laminar_pressure_ratio,
            self.atmospheric_pressure,
        ), dtype=float)
        if (not isinstance(self.cylinder, IdealDoubleActingCylinder)
                or not np.isfinite(values).all() or np.any(values <= 0)
                or self.leakage_area > self.max_area
                or not 0 < self.laminar_pressure_ratio < 1):
            raise ValueError("rectified cylinder valve settings are invalid")

    def _flow(self, pressure_drop: float, area: float,
              average_pressure: float) -> float:
        critical = ((self.atmospheric_pressure + average_pressure)
                    * (1 - self.laminar_pressure_ratio))
        return (self.discharge_coefficient * area
                * np.sqrt(2 / self.fluid_density) * pressure_drop
                / (pressure_drop * pressure_drop + critical * critical) ** .25)

    def chamber_pressures(self, rod_speed: float, high_pressure: float,
                          low_pressure: float = 0.) -> tuple[float, float]:
        """Solve A/B pressures from the two incompressible chamber balances."""
        values = np.asarray((rod_speed, high_pressure, low_pressure), dtype=float)
        if not np.isfinite(values).all() or high_pressure < low_pressure:
            raise ValueError("rod speed and ordered line pressures must be finite")

        if rod_speed > 0:
            area_a_high, area_a_low = self.leakage_area, self.max_area
            area_b_high, area_b_low = self.max_area, self.leakage_area
        else:
            area_a_high, area_a_low = self.max_area, self.leakage_area
            area_b_high, area_b_low = self.leakage_area, self.max_area

        def solve(target_flow, high_area, low_area):
            def balance(pressure):
                return (
                    self._flow(pressure - high_pressure, high_area,
                               (pressure + high_pressure) / 2)
                    + self._flow(pressure - low_pressure, low_area,
                                 (pressure + low_pressure) / 2)
                    - target_flow
                )

            span = max(high_pressure - low_pressure, 1e6)
            lower, upper = low_pressure - span, high_pressure + span
            for _ in range(32):
                if balance(lower) <= 0 and balance(upper) >= 0:
                    break
                span *= 2
                lower, upper = low_pressure - span, high_pressure + span
            else:
                raise ValueError("could not bracket cylinder chamber pressure")
            return brentq(balance, lower, upper)

        pressure_a = solve(-self.cylinder.area_a * rod_speed,
                           area_a_high, area_a_low)
        pressure_b = solve(self.cylinder.area_b * rod_speed,
                           area_b_high, area_b_low)
        return pressure_a, pressure_b

    def rod_force(self, rod_speed: float, high_pressure: float,
                  low_pressure: float = 0.) -> float:
        """Cylinder rod force in newtons from the solved chamber pressures."""
        pressure_a, pressure_b = self.chamber_pressures(
            rod_speed, high_pressure, low_pressure,
        )
        return float(self.cylinder.force(pressure_a, pressure_b))

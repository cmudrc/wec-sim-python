"""Incompressible orifice reaction used by the published OWC GBM case."""

from dataclasses import dataclass
from numbers import Real

import numpy as np


@dataclass(frozen=True)
class OrificeResponse:
    flow_rate: np.ndarray
    pressure_drop: np.ndarray
    force: np.ndarray
    absorbed_power: np.ndarray
    mach_number: np.ndarray
    compressibility_flag: np.ndarray


@dataclass(frozen=True)
class OrificePTO:
    """Evaluate the OWC piston/orifice force at a prescribed piston speed.

    ``pressure_drop`` is signed pascals, ``force`` is newtons opposing piston
    motion, and ``absorbed_power`` is positive watts. The published Simulink
    block logs pressure in kPa and the numerical power in kW. Its Mach flag
    warns when the incompressible law is outside the configured threshold;
    it does not change that law.
    """

    piston_area: float
    orifice_area: float
    discharge_coefficient: float = 0.62
    air_density: float = 1.2
    mach_threshold: float = 0.3
    sound_speed: float = 343.2

    def __post_init__(self) -> None:
        positive = (self.piston_area, self.orifice_area,
                    self.discharge_coefficient, self.air_density,
                    self.sound_speed)
        if (any(isinstance(value, bool) or not isinstance(value, Real)
                for value in (*positive, self.mach_threshold))
                or not np.isfinite((*positive, self.mach_threshold)).all()
                or any(value <= 0 for value in positive)
                or self.mach_threshold < 0):
            raise ValueError("orifice areas, discharge coefficient, air density, and sound speed must be positive; Mach threshold must be nonnegative")

    def evaluate(self, piston_speed: float | np.ndarray) -> OrificeResponse:
        speed = np.asarray(piston_speed, dtype=float)
        if not np.isfinite(speed).all():
            raise ValueError("piston speed must be finite")
        flow = self.piston_area * speed
        jet_speed = flow / self.orifice_area
        mach = jet_speed / self.sound_speed
        discharge_speed = jet_speed / self.discharge_coefficient
        pressure = 0.5 * self.air_density * discharge_speed * np.abs(discharge_speed)
        force = -self.piston_area * pressure
        return OrificeResponse(
            flow_rate=flow,
            pressure_drop=pressure,
            force=force,
            absorbed_power=np.abs(flow * pressure),
            mach_number=mach,
            compressibility_flag=np.abs(mach) > self.mach_threshold,
        )

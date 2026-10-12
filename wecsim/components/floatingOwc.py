"""Air-chamber dynamics for the published floating oscillating-water-column case."""

from dataclasses import dataclass
from numbers import Real

import numpy as np
from scipy.spatial.transform import Rotation


def _six_component_history(value, label):
    values = np.asarray(value, dtype=float)
    single = values.ndim == 1
    if single:
        values = values[None, :]
    if values.ndim != 2 or values.shape[1] != 6 or not np.isfinite(values).all():
        raise ValueError(f"{label} must contain finite six-component body states")
    return values, single


def _joint_history(value, count, label):
    values = np.asarray(value, dtype=float)
    if values.ndim == 0 and count == 1:
        values = values[None]
    if values.shape != (count,) or not np.isfinite(values).all():
        raise ValueError(f"{label} must have one finite value per body state")
    return values


@dataclass(frozen=True)
class FloatingOwcColumnJoint:
    """Rigid floater and coaxial water column joined by one axial slider.

    ``center_separation`` is the equilibrium column-center height minus the
    floater-center height. Poses contain world-center xyz and xyz Euler angles;
    velocities contain world-center velocity and world angular velocity.
    The slider's position and speed are measured in the floater's local z axis.
    """

    center_separation: float

    def __post_init__(self):
        if (isinstance(self.center_separation, bool)
                or not isinstance(self.center_separation, Real)
                or not np.isfinite(self.center_separation)
                or self.center_separation <= 0):
            raise ValueError("column center separation must be finite and positive")

    def column_pose(self, floater_pose, stroke):
        """Return the column's world-center pose from floater pose and stroke."""
        pose, single = _six_component_history(floater_pose, "floater_pose")
        slide = _joint_history(stroke, len(pose), "stroke")
        axis = Rotation.from_euler("xyz", pose[:, 3:6]).apply([0, 0, 1])
        result = pose.copy()
        result[:, :3] += (self.center_separation + slide)[:, None] * axis
        return result[0] if single else result

    def column_velocity(self, floater_pose, floater_velocity,
                        stroke, stroke_speed):
        """Return center velocity and world angular speed for the column."""
        pose, single = _six_component_history(floater_pose, "floater_pose")
        velocity, velocity_single = _six_component_history(
            floater_velocity, "floater_velocity")
        if velocity.shape != pose.shape or velocity_single != single:
            raise ValueError("floater pose and velocity histories must align")
        slide = _joint_history(stroke, len(pose), "stroke")
        slide_speed = _joint_history(stroke_speed, len(pose), "stroke_speed")
        axis = Rotation.from_euler("xyz", pose[:, 3:6]).apply([0, 0, 1])
        result = velocity.copy()
        result[:, :3] += (slide_speed[:, None] * axis
                          + (self.center_separation + slide)[:, None]
                          * np.cross(velocity[:, 3:6], axis))
        return result[0] if single else result


@dataclass(frozen=True)
class FloatingOwcChamber:
    """Compressible chamber driven by the published column heave signal.

    The published application supplies column-center world heave relative to
    its equilibrium height, and column-center world heave speed. Positive
    heave moves the water column into the chamber. ``pressure`` is gauge pressure
    in Pa, and ``turbine_speed`` is in
    rad/s. The returned force acts on the water column in the opposite direction
    from positive displacement. This is the chamber law in the published
    ``OWC/FloatingOWC`` Simulink application; body and turbine dynamics are
    separate.
    """

    area: float
    initial_volume: float
    gamma: float
    ambient_pressure: float
    ambient_density: float
    turbine_diameter: float
    turbine_kappa: float

    def __post_init__(self) -> None:
        values = tuple(self.__dict__.values())
        if (any(isinstance(value, bool) or not isinstance(value, Real)
                for value in values)
                or not np.isfinite(values).all()
                or any(value <= 0 for value in values)):
            raise ValueError("floating OWC chamber parameters must be finite and positive")

    def volume(self, displacement: float | np.ndarray) -> np.ndarray:
        volume = self.initial_volume - self.area * np.asarray(displacement, dtype=float)
        if not np.isfinite(volume).all() or np.any(volume <= 0):
            raise ValueError("floating OWC chamber volume must remain positive")
        return volume

    def pressure_derivative(
        self,
        pressure: float | np.ndarray,
        displacement: float | np.ndarray,
        speed: float | np.ndarray,
        turbine_speed: float | np.ndarray,
    ) -> np.ndarray:
        """Return the gauge-pressure derivative in Pa/s."""
        gauge = np.asarray(pressure, dtype=float)
        speed = np.asarray(speed, dtype=float)
        omega = np.asarray(turbine_speed, dtype=float)
        absolute = self.ambient_pressure + gauge
        if (not np.isfinite(gauge).all() or not np.isfinite(speed).all()
                or not np.isfinite(omega).all() or np.any(absolute <= 0)
                or np.any(omega == 0)):
            raise ValueError("pressure, water-column speed, and nonzero turbine speed must be physical")
        density = self.ambient_density * (
            absolute / self.ambient_pressure) ** (1 / self.gamma)
        turbine_volume_flow = (
            gauge * self.turbine_diameter / (self.turbine_kappa * omega * density)
        )
        return self.gamma * absolute / self.volume(displacement) * (
            self.area * speed - turbine_volume_flow
        )

    def force_on_column(self, pressure: float | np.ndarray) -> np.ndarray:
        """Return the chamber's axial force on the water column in N."""
        gauge = np.asarray(pressure, dtype=float)
        if not np.isfinite(gauge).all():
            raise ValueError("pressure must be finite")
        return -self.area * gauge


def _published_wells_curves() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Sample the pinned application's fitted Wells curves on its 5000-point grid."""
    psi = np.linspace(-0.3, 0.3, 5000)
    magnitude = np.abs(psi)
    efficiency = np.zeros_like(psi)
    linear_segments = (
        (0.008, 0.009, 0.0, 0.1),
        (0.009, 0.0128, 0.1, 0.2),
        (0.0128, 0.0199, 0.2, 0.3),
        (0.0199, 0.0404, 0.3, 0.5),
        (0.0878, 0.1141, 0.5, 0.3),
        (0.1141, 0.1192, 0.3, 0.273),
        (0.2122, 0.3, 0.0016, 0.0),
    )
    for low_psi, high_psi, low_eta, high_eta in linear_segments:
        selected = (magnitude >= low_psi) & (magnitude <= high_psi)
        efficiency[selected] = low_eta + (high_eta - low_eta) * (
            magnitude[selected] - low_psi) / (high_psi - low_psi)
    coefficients = np.polyfit(
        [0.0404, 0.0634, 0.0878], [0.5, 0.595, 0.5], 2)
    selected = (magnitude > 0.0404) & (magnitude < 0.0878)
    efficiency[selected] = np.polyval(coefficients, magnitude[selected])
    selected = (magnitude >= 0.1192) & (magnitude <= 0.2122)
    efficiency[selected] = 200.6 * np.exp(-55.35 * magnitude[selected])

    flow = np.where(
        magnitude <= 0.1, 0.775 * magnitude,
        -2.503 * magnitude**2 + 1.608 * magnitude - 0.05664 - 0.00167,
    )
    torque = efficiency * magnitude * flow
    selected = (magnitude > 0.1) & (magnitude <= 0.114)
    torque[selected] = (0.3156 + (0.3170 - 0.3156) *
                        (magnitude[selected] - 0.1) / 0.014) / 100
    return psi, efficiency, flow, torque


_PSI_GRID, _EFFICIENCY, _FLOW, _TORQUE = _published_wells_curves()


@dataclass(frozen=True)
class FloatingOwcTurbineResponse:
    control_torque: np.ndarray
    turbine_torque: np.ndarray
    load_power: np.ndarray
    pneumatic_power: np.ndarray
    efficiency: np.ndarray
    speed_derivative: np.ndarray


@dataclass(frozen=True)
class FloatingOwcTurbine:
    """Published floating-OWC Wells turbine, load control, and rotor inertia.

    ``evaluate`` accepts chamber gauge pressure in Pa and rotor speed in
    rad/s. ``load_power`` is the source's logged ``P_turb = u * speed`` in W;
    ``pneumatic_power`` is its pressure-flow power in W. The performance
    curves are the fixed fits from the published application.
    """

    diameter: float = 0.75
    inertia: float = 3.06
    ambient_pressure: float = 101325.0
    ambient_density: float = 1.25
    gamma: float = 1.4
    max_speed: float = 350.0
    control_coefficient: float = 2e-4
    control_exponent: float = 3.0
    max_control_torque: float = 216.5

    def __post_init__(self) -> None:
        values = tuple(self.__dict__.values())
        if (any(isinstance(value, bool) or not isinstance(value, Real)
                for value in values)
                or not np.isfinite(values).all()
                or any(value <= 0 for value in values)
                or self.control_exponent < 1):
            raise ValueError("floating OWC turbine parameters must be finite and positive; control exponent must be at least one")

    def evaluate(
        self, pressure: float | np.ndarray, speed: float | np.ndarray,
    ) -> FloatingOwcTurbineResponse:
        gauge = np.asarray(pressure, dtype=float)
        rotor_speed = np.asarray(speed, dtype=float)
        absolute = self.ambient_pressure + gauge
        if (not np.isfinite(gauge).all() or not np.isfinite(rotor_speed).all()
                or np.any(absolute <= 0) or np.any(rotor_speed <= 0)):
            raise ValueError("chamber absolute pressure and rotor speed must be positive")

        chamber_density = self.ambient_density * (
            absolute / self.ambient_pressure) ** (1 / self.gamma)
        density = np.maximum(self.ambient_density, chamber_density)
        pressure_coefficient = gauge / (
            density * rotor_speed**2 * self.diameter**2)
        safe_psi = np.where(rotor_speed > self.max_speed,
                            0.0, pressure_coefficient)
        clipped = np.clip(safe_psi, _PSI_GRID[0], _PSI_GRID[-1])
        flow = np.interp(clipped, _PSI_GRID, _FLOW) * np.sign(pressure_coefficient)
        torque_coefficient = np.interp(clipped, _PSI_GRID, _TORQUE)
        efficiency = np.interp(clipped, _PSI_GRID, _EFFICIENCY)

        turbine_torque = (self.diameter**5 * density * torque_coefficient *
                          rotor_speed**2)
        control_torque = np.minimum(
            self.control_coefficient * rotor_speed**(self.control_exponent - 1),
            self.max_control_torque)
        mass_flow = flow * density * rotor_speed * self.diameter**3
        return FloatingOwcTurbineResponse(
            control_torque=control_torque,
            turbine_torque=turbine_torque,
            load_power=control_torque * rotor_speed,
            pneumatic_power=gauge * mass_flow / chamber_density,
            efficiency=efficiency,
            speed_derivative=(turbine_torque - control_torque) / self.inertia,
        )

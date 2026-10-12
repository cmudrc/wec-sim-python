"""Discrete PTO control laws for published WEC-Sim application cases."""

from dataclasses import dataclass
import math


@dataclass(frozen=True)
class DeclutchingState:
    """Controller memory before the next velocity sample."""

    elapsed_off: float = 0.0
    elapsed_on: float = 0.0
    previous_velocity: float = 0.0


@dataclass(frozen=True)
class DeclutchingControl:
    """Zero PTO force for a timed interval after a velocity sign change.

    This is the discrete law in the pinned Sphere Declutching application.
    It samples the PTO velocity once per simulation step. ``gain`` is the
    proportional damping applied while engaged.
    """

    gain: float
    declutch_time: float
    minimum_on_time: float = 0.2

    def __post_init__(self):
        if not all(math.isfinite(value) for value in (
                self.gain, self.declutch_time, self.minimum_on_time)):
            raise ValueError("declutching settings must be finite")
        if self.gain <= 0 or self.declutch_time <= 0 or self.minimum_on_time < 0:
            raise ValueError("declutching settings are outside their supported range")

    def initial_state(self) -> DeclutchingState:
        return DeclutchingState()

    def sample(self, velocity: float, state: DeclutchingState,
               dt: float) -> tuple[float, DeclutchingState]:
        """Return signed force and next memory for one major simulation step."""
        if not isinstance(state, DeclutchingState):
            raise TypeError("state must be DeclutchingState")
        if not (math.isfinite(velocity) and math.isfinite(dt) and dt > 0):
            raise ValueError("velocity must be finite and dt positive")
        reversal = (self._sign(velocity) != self._sign(state.previous_velocity)
                    and state.elapsed_off == 0
                    and state.elapsed_on > self.minimum_on_time)
        interval_active = 0 < state.elapsed_off < self.declutch_time
        if reversal or interval_active:
            return 0.0, DeclutchingState(
                state.elapsed_off + dt, 0.0, velocity,
            )
        return -self.gain * velocity, DeclutchingState(
            0.0, state.elapsed_on + dt, velocity,
        )

    @staticmethod
    def _sign(value: float) -> int:
        return (value > 0) - (value < 0)


@dataclass(frozen=True)
class LatchingState:
    """Elapsed latching and normal intervals before the next velocity sample."""

    elapsed_latched: float = 0.0
    elapsed_normal: float = 0.0
    previous_velocity: float = 0.0


@dataclass(frozen=True)
class LatchingControl:
    """Apply strong damping after a velocity reversal for a timed interval.

    The pinned Sphere Latching application implements a finite, dissipative
    damping force during the latch; it does not impose a rigid displacement
    constraint. Controller memory advances once per major simulation step.
    """

    gain: float
    latch_damping: float
    latch_time: float
    minimum_normal_time: float = 0.2

    def __post_init__(self):
        values = (self.gain, self.latch_damping, self.latch_time,
                  self.minimum_normal_time)
        if not all(math.isfinite(value) for value in values):
            raise ValueError("latching settings must be finite")
        if (self.gain <= 0 or self.latch_damping <= 0
                or self.latch_time <= 0 or self.minimum_normal_time < 0):
            raise ValueError("latching settings are outside their supported range")

    def initial_state(self) -> LatchingState:
        return LatchingState()

    def sample(self, velocity: float, state: LatchingState,
               dt: float) -> tuple[float, LatchingState]:
        """Return signed force and next memory for one major simulation step."""
        if not isinstance(state, LatchingState):
            raise TypeError("state must be LatchingState")
        if not (math.isfinite(velocity) and math.isfinite(dt) and dt > 0):
            raise ValueError("velocity must be finite and dt positive")
        reversal = (DeclutchingControl._sign(velocity)
                    != DeclutchingControl._sign(state.previous_velocity)
                    and state.elapsed_latched == 0
                    and state.elapsed_normal > self.minimum_normal_time)
        interval_active = 0 < state.elapsed_latched < self.latch_time
        if reversal or interval_active:
            return -self.latch_damping * velocity, LatchingState(
                state.elapsed_latched + dt, 0.0, velocity,
            )
        return -self.gain * velocity, LatchingState(
            0.0, state.elapsed_normal + dt, velocity,
        )

"""WEC builder, configuration records, and simulation results.

``from wecsim import WEC, RegularWave`` is the concise public import. This
package also groups the Python-facing modeling types for explicit imports.
"""

from .api import WEC
from .configuration import (
    Body,
    BodyPoint,
    Coordinate,
    Current,
    FullDirectionalSpectrumWave,
    HydroState,
    ImportedElevationWave,
    ImportedSpectrumWave,
    JONSWAPWave,
    LinearPTO,
    Motion,
    MultiPMWave,
    NoWave,
    PMWave,
    RegularCICWave,
    RegularWave,
    RotationalPTO,
    SimpleDirectDrive,
    VariableHydro,
    WorldPoint,
)
from .results import (
    DirectDriveHistory,
    FlexibleModeHistory,
    LinearGeneratorHistory,
    MotionHistory,
    PTOHistory,
    WECResult,
)

__all__ = [
    "Body", "BodyPoint", "Coordinate", "Current", "DirectDriveHistory",
    "FlexibleModeHistory", "FullDirectionalSpectrumWave", "HydroState",
    "ImportedElevationWave", "ImportedSpectrumWave", "JONSWAPWave",
    "LinearGeneratorHistory", "LinearPTO", "Motion", "MotionHistory",
    "MultiPMWave", "NoWave", "PMWave", "PTOHistory", "RegularCICWave",
    "RegularWave", "RotationalPTO", "SimpleDirectDrive", "VariableHydro",
    "WEC", "WECResult", "WorldPoint",
]

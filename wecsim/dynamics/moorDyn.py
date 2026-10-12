"""Optional binding to the legacy MoorDyn library used by MATLAB WEC-Sim.

The native library and its input file are supplied by the caller. MoorDyn v2's
legacy ``MoorDynInit``/``MoorDynStep`` API has one process-wide simulation.
"""

import ctypes
import os
from pathlib import Path
from threading import Lock

import numpy as np


_active_session = Lock()
_Six = ctypes.c_double * 6


def _six(values, name):
    vector = np.asarray(values, dtype=float)
    if vector.shape != (6,) or not np.isfinite(vector).all():
        raise ValueError(f"{name} must be a finite six-component vector")
    return _Six(*vector)


class MoorDyn:
    """Step one externally supplied MoorDyn v2 mooring model from Python.

    Start with a six-component world pose and velocity, then call ``step``
    once per coupling interval. The returned six-vector is the force and
    moment exerted by the lines on the coupled body. Use ``with`` to ensure
    the native model closes. The input file controls line geometry and writes
    MoorDyn's diagnostic output beside itself.
    """

    def __init__(self, library_file: str | Path, input_file: str | Path):
        self.library_file = Path(library_file).expanduser().resolve(strict=True)
        self.input_file = Path(input_file).expanduser().resolve(strict=True)
        self._library = None
        self._started = False

    def start(self, position, velocity):
        """Initialize the process-wide native model and return this session."""
        if self._started:
            raise RuntimeError("MoorDyn is already started")
        initial_position = _six(position, "position")
        initial_velocity = _six(velocity, "velocity")
        if not _active_session.acquire(blocking=False):
            raise RuntimeError("another MoorDyn session is active")
        try:
            library = ctypes.CDLL(str(self.library_file))
            library.SetDisableConsole.argtypes = [ctypes.c_int]
            library.SetDisableConsole.restype = None
            library.MoorDynInit.argtypes = [ctypes.POINTER(ctypes.c_double),
                                           ctypes.POINTER(ctypes.c_double),
                                           ctypes.c_char_p]
            library.MoorDynInit.restype = ctypes.c_int
            library.MoorDynStep.argtypes = [ctypes.POINTER(ctypes.c_double),
                                           ctypes.POINTER(ctypes.c_double),
                                           ctypes.POINTER(ctypes.c_double),
                                           ctypes.POINTER(ctypes.c_double),
                                           ctypes.POINTER(ctypes.c_double)]
            library.MoorDynStep.restype = ctypes.c_int
            library.MoorDynClose.argtypes = []
            library.MoorDynClose.restype = ctypes.c_int
            library.SetDisableConsole(1)
            code = library.MoorDynInit(initial_position, initial_velocity,
                                       os.fsencode(self.input_file))
            if code != 0:
                library.MoorDynClose()
                raise RuntimeError(f"MoorDynInit failed with code {code}")
            self._library = library
            self._started = True
            return self
        except BaseException:
            _active_session.release()
            raise

    def step(self, position, velocity, time: float, dt: float) -> np.ndarray:
        """Advance the mooring from ``time`` by ``dt`` and return its load."""
        if not self._started:
            raise RuntimeError("MoorDyn must be started before stepping")
        pose = _six(position, "position")
        speed = _six(velocity, "velocity")
        if not np.isfinite([time, dt]).all() or dt <= 0:
            raise ValueError("time and positive dt must be finite")
        native_time = ctypes.c_double(time)
        native_dt = ctypes.c_double(dt)
        force = _Six(*([0.0] * 6))
        code = self._library.MoorDynStep(pose, speed, force,
                                         ctypes.byref(native_time),
                                         ctypes.byref(native_dt))
        if code != 0:
            raise RuntimeError(f"MoorDynStep failed with code {code}")
        return np.asarray(force[:], dtype=float)

    def close(self):
        """Close the native model; repeated calls are harmless."""
        if self._started:
            self._started = False
            try:
                code = self._library.MoorDynClose()
            finally:
                self._library = None
                _active_session.release()
            if code != 0:
                raise RuntimeError(f"MoorDynClose failed with code {code}")

    def __enter__(self):
        if not self._started:
            raise RuntimeError("call start before entering the MoorDyn context")
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

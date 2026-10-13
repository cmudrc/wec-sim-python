"""Check the legacy MoorDyn call boundary without requiring a native install."""

import ctypes
import os
from xml.etree import ElementTree

import numpy as np
import pytest

from wecsim.dynamics.moorDyn import MoorDyn
from wecsim.visualization.paraviewClass import ParaviewClass


class _Function:
    def __init__(self, callback):
        self.callback = callback

    def __call__(self, *args):
        return self.callback(*args)


class _NativeMoorDyn:
    def __init__(self):
        self.calls = []
        self.load = np.arange(6, dtype=float) + 1
        self.init_code = 0
        self.step_code = 0
        self.SetDisableConsole = _Function(self._disable_console)
        self.MoorDynInit = _Function(self._init)
        self.MoorDynStep = _Function(self._step)
        self.MoorDynClose = _Function(self._close)

    def _disable_console(self, value):
        self.calls.append(("console", value))

    def _init(self, pose, speed, input_file):
        self.calls.append(("init", np.array(pose[:]), np.array(speed[:]),
                           os.fsdecode(input_file)))
        return self.init_code

    def _step(self, pose, speed, load, time, dt):
        self.calls.append(("step", np.array(pose[:]), np.array(speed[:]),
                           ctypes.cast(time, ctypes.POINTER(ctypes.c_double)).contents.value,
                           ctypes.cast(dt, ctypes.POINTER(ctypes.c_double)).contents.value))
        for index, value in enumerate(self.load):
            load[index] = value
        return self.step_code

    def _close(self):
        self.calls.append(("close",))
        return 0


def test_native_six_vector_contract_and_exclusive_session(tmp_path, monkeypatch):
    library_file = tmp_path / "libmoordyn.so"
    input_file = tmp_path / "lines.txt"
    library_file.touch()
    input_file.touch()
    native = _NativeMoorDyn()
    monkeypatch.setattr(ctypes, "CDLL", lambda path: native)
    first = MoorDyn(library_file, input_file)
    second = MoorDyn(library_file, input_file)
    pose = np.arange(6, dtype=float)
    speed = pose / 10
    with first.start(pose, speed) as session:
        with pytest.raises(RuntimeError, match="another MoorDyn session"):
            second.start(pose, speed)
        np.testing.assert_array_equal(session.step(pose + 1, speed, 0.25, 0.01),
                                      native.load)
        assert native.calls[1][0] == "init"
        np.testing.assert_array_equal(native.calls[1][1], pose)
        np.testing.assert_array_equal(native.calls[1][2], speed)
        assert native.calls[1][3] == str(input_file)
        assert native.calls[2][0] == "step"
        np.testing.assert_array_equal(native.calls[2][1], pose + 1)
        np.testing.assert_array_equal(native.calls[2][2], speed)
        assert native.calls[2][3:] == (0.25, 0.01)
        with pytest.raises(ValueError, match="finite six-component"):
            session.step(pose, np.full(6, np.nan), 0.25, 0.01)
    assert native.calls[-1] == ("close",)
    with second.start(pose, speed):
        pass
    assert native.calls[-1] == ("close",)


def test_native_failures_are_reported_without_leaking_session(tmp_path, monkeypatch):
    library_file = tmp_path / "libmoordyn.so"
    input_file = tmp_path / "lines.txt"
    library_file.touch()
    input_file.touch()
    native = _NativeMoorDyn()
    monkeypatch.setattr(ctypes, "CDLL", lambda path: native)
    mooring = MoorDyn(library_file, input_file)
    zero = np.zeros(6)

    native.init_code = 7
    with pytest.raises(RuntimeError, match="MoorDynInit failed with code 7"):
        mooring.start(zero, zero)
    assert native.calls[-1] == ("close",)

    native.init_code = 0
    with mooring.start(zero, zero):
        native.step_code = 9
        with pytest.raises(RuntimeError, match="MoorDynStep failed with code 9"):
            mooring.step(zero, zero, 0, 0.01)
        native.step_code = 0
        native.load[4] = np.nan
        with pytest.raises(RuntimeError, match="nonfinite connection load"):
            mooring.step(zero, zero, 0, 0.01)
    assert native.calls[-1] == ("close",)


def test_native_line_history_can_feed_mooring_vtp(tmp_path):
    """Use MoorDyn's line output directly to render its geometry and tension."""
    library_file = tmp_path / "libmoordyn.so"
    input_file = tmp_path / "lines.txt"
    library_file.touch()
    input_file.touch()
    model = MoorDyn(library_file, input_file)
    output = tmp_path / "lines_Line1.out"
    with pytest.raises(FileNotFoundError, match="LineOutputs"):
        model.line_history(1)
    output.write_text(
        "Time Node0px Node0py Node0pz Node1px Node1py Node1pz Seg1Te\n"
        "0 0 0 -10 1 0 -5 20\n"
        "1 0 0 -10 2 0 -5 30\n",
        encoding="utf-8",
    )
    times, nodes, tensions = model.line_history(1)
    np.testing.assert_array_equal(times, [0, 1])
    assert nodes.shape == (2, 2, 3)
    assert tensions.shape == (2, 1)
    paths = ParaviewClass(None).write_paraview_vtp_mooring(
        times, [0, 0.5, 1], tmp_path / "vtp",
        lines=[nodes], tensions=[tensions],
    )
    middle = ElementTree.parse(paths[1]).find("./PolyData/Piece")
    points = np.fromstring(middle.findtext("./Points/DataArray"), sep=" ")
    load = np.fromstring(
        middle.findtext("./CellData/DataArray[@Name='Segment Tension']"),
        sep=" ",
    )
    np.testing.assert_allclose(points.reshape(2, 3),
                               [[0, 0, -10], [1.5, 0, -5]])
    np.testing.assert_allclose(load, [25])

"""Reject passive-yaw configurations whose directional physics are absent."""

from pathlib import Path

import pytest

from wecsim import PMWave, RegularWave, WEC


HYDRO = (Path(__file__).parent / "test_objects/test_bodyclass/testData"
         / "hydroData/oswec.h5")


def test_passive_yaw_rejects_non_yaw_motion():
    wec = WEC("Unsupported passive yaw")
    flap = wec.body("flap", HYDRO, mass=12700,
                    inertia=(1.85e6,) * 3, passive_yaw=True)
    wec.coordinate("heave", flap.move("heave"))
    with pytest.raises(ValueError, match="one pure-yaw body"):
        wec.run(RegularWave(2.5, 8, 10), dt=0.1,
                end_time=0.1, ramp_time=100)


def test_passive_yaw_rejects_incomplete_bem_heading_range():
    # This historical OSWEC fixture contains headings only from 0 to 90°.
    wec = WEC("Incomplete passive yaw")
    flap = wec.body("flap", HYDRO, mass=12700,
                    inertia=(1.85e6,) * 3, passive_yaw=True)
    wec.coordinate("yaw", flap.move("yaw"))
    with pytest.raises(ValueError, match="full-circle BEM headings"):
        wec.run(RegularWave(2.5, 8, 10), dt=0.1,
                end_time=0.1, ramp_time=100)


@pytest.mark.parametrize("threshold", [-1, float("nan"), True, "1"])
def test_passive_yaw_rejects_invalid_threshold(threshold):
    with pytest.raises(ValueError, match="passive_yaw_threshold"):
        WEC().body("flap", HYDRO, passive_yaw=True,
                   passive_yaw_threshold=threshold)


def test_heading_threshold_requires_passive_yaw_and_pm_wave():
    with pytest.raises(ValueError, match="passive_yaw_threshold"):
        WEC().body("flap", HYDRO, passive_yaw_threshold=1)
    wec = WEC("PM-only heading threshold")
    flap = wec.body("flap", HYDRO, mass=12700,
                    inertia=(1.85e6,) * 3, passive_yaw=True,
                    passive_yaw_threshold=1)
    wec.coordinate("yaw", flap.move("yaw"))
    with pytest.raises(ValueError, match="positive passive_yaw_threshold needs PM"):
        wec.run(RegularWave(2.5, 8, 10), dt=0.1,
                end_time=0.1, ramp_time=100)


def test_heading_bank_accepts_pm_wave_before_hydro_validation():
    wec = WEC("PM heading bank")
    flap = wec.body("flap", HYDRO, mass=12700,
                    inertia=(1.85e6,) * 3, passive_yaw=True,
                    yaw_heading_bank=[-10, 0, 10])
    wec.coordinate("yaw", flap.move("yaw"))
    # This small historical fixture is intentionally missing full-circle
    # headings. The PM bank must get past configuration validation first.
    with pytest.raises(ValueError, match="full-circle BEM headings"):
        wec.run(PMWave(2.5, 8, direction=10), dt=0.1,
                end_time=0.1, ramp_time=0)


def test_source_ode4_rejects_devices_without_held_passive_yaw():
    wec = WEC("Unsupported source solver")
    body = wec.body("flap", HYDRO, mass=12700,
                    inertia=(1.85e6,) * 3)
    wec.coordinate("yaw", body.move("yaw"))
    with pytest.raises(ValueError, match="PM held yaw or a heading bank"):
        wec.run(PMWave(2.5, 8, direction=10), dt=.1,
                end_time=.1, ramp_time=0, solver="wecsim_ode4")

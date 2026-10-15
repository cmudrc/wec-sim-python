"""Checks against MATLAB-generated fixtures and the WEC-Sim wave equations.

The older object tests import copies of the implementation from their own
directories. These tests deliberately import the production WaveClass.
"""

from pathlib import Path
import sys
from xml.etree import ElementTree

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]

from wecsim.environment.waveClass import WaveClass  # noqa: E402
from wecsim.visualization.paraviewClass import ParaviewClass  # noqa: E402


def test_regular_wave_matches_matlab_fixture():
    wave = WaveClass("regular")
    wave.T = 8.0
    wave.H = 2.5
    wave.waveDir = [0]
    wave.wavegauge1loc = [5, 5]
    wave.wavegauge2loc = [10, 0]
    wave.wavegauge3loc = [0, -10]
    wave.waveSetup(
        [5.19999512307279, 0.0199999977946844],
        "infinite", 100.0, 0.1, 2000, 9.81, 1000.0, 200.0,
    )

    fixture = ROOT / "tests/test_objects/test_waveclass/testData/regular_1_test"
    for attribute, filename in (
        ("waveAmpTime", "waveAmpTime.txt"),
        ("waveAmpTime1", "waveAmpTime1.txt"),
        ("waveAmpTime2", "waveAmpTime2.txt"),
        ("waveAmpTime3", "waveAmpTime3.txt"),
    ):
        expected = np.loadtxt(fixture / filename).T
        np.testing.assert_allclose(getattr(wave, attribute), expected, rtol=0, atol=1e-12)


def test_zero_ramp_directional_gauges_match_current_matlab():
    """Zero ramp and finite depth keep the spatial phase at each marker."""
    wave = WaveClass("regular")
    wave.T = 8
    wave.H = 2.5
    wave.waveDir = [30]
    wave.wavegauge1loc = [4, 0]
    wave.wavegauge2loc = [0, 5]
    wave.wavegauge3loc = [-3, 2]
    wave.waveSetup([0.4, 2], 25, 0, 0.1, 20, 9.81, 1000, 2)

    source = np.loadtxt(
        ROOT / "tests/fixtures/matlab/zero_ramp_regular_gauges.csv",
        delimiter=",",
    )
    assert source.shape == (21, 5)
    np.testing.assert_allclose(source[:, 0], np.arange(21) * 0.1,
                               rtol=0, atol=1e-14)
    actual = np.column_stack((
        wave.waveAmpTime[0], wave.waveAmpTime[1],
        wave.waveAmpTime1[1], wave.waveAmpTime2[1], wave.waveAmpTime3[1],
    ))
    np.testing.assert_allclose(actual, source, rtol=0, atol=1e-12)


def test_finite_depth_regular_power_follows_matlab_group_velocity():
    wave = WaveClass("regular")
    wave.w = 2 * np.pi / 8
    wave.A = 1.25
    wave.deepWaterWave = 0
    wave.waterDepth = 25.0
    wave.waveNumber(9.81)
    wave.wavePowerReg(9.81, 1000.0)

    k = wave.k
    np.testing.assert_allclose(wave.w**2, 9.81 * k * np.tanh(k * 25.0), rtol=1e-13)
    phase_speed = wave.w / k
    group_speed = phase_speed / 2 * (1 + 2 * k * 25.0 / np.sinh(2 * k * 25.0))
    expected_power = 0.5 * 1000.0 * 9.81 * wave.A**2 * group_speed
    np.testing.assert_allclose(wave.Pw, expected_power, rtol=1e-13)


def test_historical_equal_energy_matlab_reference():
    wave = WaveClass("irregular")
    wave.w = np.array([
        0.525743641856087, 0.541511547017433, 0.551933697209493,
        0.559973049643924, 0.566624163384781, 0.572373957973840,
        0.577471073177114, 0.582081268838611, 0.586297784870588,
        0.590203501195047, 0.593850217763243, 0.597279374536177,
        0.600522051484601, 0.603609328579267, 0.606551565810425,
    ])
    wave.dw = np.array([
        0.505743644061402, 0.0157679051613465, 0.0104221501920596,
        0.00803935243443155, 0.00665111374085714, 0.00574979458905867,
        0.00509711520327372, 0.00461019566149745, 0.00421651603197637,
        0.00390571632445969, 0.00364671656819582, 0.00342915677293409,
        0.00324267694842406, 0.00308727709466572, 0.00294223723115805,
    ])
    wave.T = 8
    wave.H = 2.5
    wave.spectrumType = "BS"
    wave.freqDisc = "EqualEnergy"
    wave.deepWaterWave = 1
    wave.numFreq = 4
    wave.irregWaveSpectrum(9.81, 1000.0)

    np.testing.assert_allclose(
        wave.w, [0.566624163384781, 0.582081268838611,
                 0.590203501195047, 0.600522051484601], rtol=1e-13,
    )
    np.testing.assert_allclose(
        wave.S, [0.126744401973820, 0.177305360376672,
                 0.206784788616044, 0.246458761526240], rtol=1e-8,
    )
    np.testing.assert_allclose(wave.Pw, 16836.8900561635, rtol=1e-13)


def test_irregular_wave_setup_uses_height_and_local_seed():
    np.random.seed(731)
    before = np.random.get_state()
    waves = []
    for height in (2.5, 5.0, 2.5):
        wave = WaveClass("irregular")
        wave.T = 8
        wave.H = height
        wave.spectrumType = "PM"
        wave.freqDisc = "Traditional"
        wave.numFreq = 64
        wave.phaseSeed = 1
        wave.waveSetup([0.4, 2.0], "infinite", 1, 0.1, 20, 9.81, 1000, 2)
        waves.append(wave)
    after = np.random.get_state()
    assert before[0] == after[0]
    np.testing.assert_array_equal(before[1], after[1])
    assert before[2:] == after[2:]
    np.testing.assert_allclose(waves[1].S, 4 * waves[0].S, rtol=1e-14)
    np.testing.assert_array_equal(waves[0].phase, waves[2].phase)
    assert np.isfinite(waves[0].waveAmpTime[1]).all()


def test_irregular_wave_setup_initializes_traditional_default_count():
    wave = WaveClass("irregular")
    wave.T = 8
    wave.H = 2.5
    wave.spectrumType = "PM"
    wave.freqDisc = "Traditional"
    wave.phaseSeed = 1
    wave.waveSetup([0.4, 2.0], "infinite", 1, 0.1, 20, 9.81, 1000, 2)
    assert wave.numFreq == 1000
    assert wave.S.shape == (1000,)


def test_irregular_phase_replay_rejects_missing_directions():
    wave = WaveClass("irregular")
    wave.T = 8
    wave.H = 2.5
    wave.spectrumType = "PM"
    wave.freqDisc = "Traditional"
    wave.numFreq = 16
    wave.waveDir = [0, 30]
    wave.waveSpread = [0.5, 0.5]
    wave.phaseData = np.zeros((16, 1))
    with pytest.raises(ValueError, match=r"shape \(frequency, direction\)"):
        wave.waveSetup([0.4, 2.0], "infinite", 1, 0.1, 20, 9.81, 1000, 2)


@pytest.mark.parametrize("ramp_time", [0.0, 1.0])
def test_directional_irregular_gauges_follow_spatial_phase_and_ramp(ramp_time):
    """Each direction shifts a gauge by its projected travel distance."""
    wave = WaveClass("irregular")
    wave.w = np.array([1.0, 1.8])
    wave.k = np.array([0.15, 0.3])
    wave.A = np.array([0.8, 1.2])
    wave.dw = np.array([0.2, 0.3])
    wave.waveDir = [0, 90]
    wave.waveSpread = [0.7, 0.3]
    wave.phase = np.array([[0.1, 0.3], [-0.2, 0.4]])
    gauges = [(0, 0), (4, 0), (0, 3), (-1, 2)]
    wave.wavegauge1loc, wave.wavegauge2loc, wave.wavegauge3loc = gauges[1:]

    wave.waveElevIrreg(ramp_time, 0.25, 8, wave.dw)

    times = np.arange(9) * 0.25
    envelope = (np.ones_like(times) if ramp_time == 0 else
                np.sin(np.pi * np.minimum(times, ramp_time) / (2 * ramp_time))**2)
    for point, history in zip(gauges, (
            wave.waveAmpTime, wave.waveAmpTime1,
            wave.waveAmpTime2, wave.waveAmpTime3)):
        expected = np.zeros_like(times)
        for heading, spread, phases in zip(
                wave.waveDir, wave.waveSpread, wave.phase):
            projected_distance = (point[0] * np.cos(np.deg2rad(heading))
                                  + point[1] * np.sin(np.deg2rad(heading)))
            amplitudes = np.sqrt(wave.A * wave.dw * spread)
            expected += np.sum(amplitudes[:, None] * np.cos(
                wave.w[:, None] * times
                - wave.k[:, None] * projected_distance
                + phases[:, None],
            ), axis=0)
        np.testing.assert_array_equal(history[0], times)
        np.testing.assert_allclose(history[1], expected * envelope,
                                   rtol=0, atol=1e-14)


def test_wave_surface_grid_matches_matlab_regular_and_irregular_equations():
    X = np.array([[0.0, 2.0], [0.0, 2.0]])
    Y = np.array([[0.0, 0.0], [3.0, 3.0]])

    no_wave = WaveClass("noWave")
    np.testing.assert_array_equal(
        ParaviewClass(no_wave).waveElevationGrid(2.0, X, Y), np.zeros_like(X)
    )

    regular = WaveClass("regular")
    regular.A = 1.25
    regular.k = 0.3
    regular.w = 0.8
    regular.waveDir = [90]
    expected = 1.25 * np.cos(-0.3 * Y + 0.8 * 2.0)
    np.testing.assert_allclose(
        ParaviewClass(regular).waveElevationGrid(2.0, X, Y), expected, atol=1e-15
    )

    irregular = WaveClass("irregular")
    irregular.A = np.array([2.0, 1.0])
    irregular.dw = np.array([0.1, 0.2])
    irregular.w = np.array([0.8, 1.2])
    irregular.k = np.array([0.3, 0.5])
    irregular.waveDir = [0, 90]
    irregular.waveSpread = [0.6, 0.4]
    irregular.phase = np.array([[0.0, 0.1], [0.2, 0.3]])
    expected = np.zeros_like(X)
    for direction_index, position in enumerate((X, Y)):
        for frequency_index in range(2):
            expected += np.sqrt(
                irregular.A[frequency_index] * irregular.dw[frequency_index]
                * irregular.waveSpread[direction_index]
            ) * np.cos(
                -irregular.k[frequency_index] * position
                + irregular.w[frequency_index] * 2.0
                + irregular.phase[direction_index, frequency_index]
            )
    np.testing.assert_allclose(
        ParaviewClass(irregular).waveElevationGrid(2.0, X, Y), expected, atol=1e-15
    )


def test_paraview_wave_vtp_matches_published_grid_layout(tmp_path):
    wave = WaveClass("regular")
    wave.A = 1.25
    wave.k = 0.3
    wave.w = 0.8
    wave.waveDir = [90]
    wave.waterDepth = 30
    paths = ParaviewClass(wave).write_paraview_vtp_wave(
        [2.0, 3.0], tmp_path, domain_size=6,
        num_points_x=3, num_points_y=2,
    )
    assert [path.name for path in paths] == ["waves_1.vtp", "waves_2.vtp"]
    assert (tmp_path / "ground.txt").read_text() == "6\n30\n0\n"
    X, Y = np.meshgrid([-6.0, 0.0, 6.0], [-6.0, 6.0])
    for time, path in zip((2.0, 3.0), paths):
        root = ElementTree.parse(path).getroot()
        assert root.attrib == {"type": "PolyData", "version": "0.1"}
        piece = root.find("./PolyData/Piece")
        assert piece.attrib == {"NumberOfPoints": "6", "NumberOfPolys": "2"}
        points = np.fromstring(piece.findtext("./Points/DataArray"), sep=" ").reshape(-1, 3)
        expected = np.column_stack((
            X.ravel(), Y.ravel(),
            (1.25 * np.cos(-0.3 * Y + 0.8 * time)).ravel(),
        ))
        np.testing.assert_allclose(points, expected, rtol=0, atol=5.1e-6)
        connectivity = np.fromstring(
            piece.findtext("./Polys/DataArray[@Name='connectivity']"),
            sep=" ", dtype=int,
        ).reshape(-1, 4)
        np.testing.assert_array_equal(connectivity, [[0, 1, 4, 3], [1, 2, 5, 4]])
        offsets = np.fromstring(
            piece.findtext("./Polys/DataArray[@Name='offsets']"),
            sep=" ", dtype=int,
        )
        np.testing.assert_array_equal(offsets, [4, 8])
    with pytest.raises(ValueError, match="finite and increasing"):
        ParaviewClass(wave).write_paraview_vtp_wave(
            [2, 2], tmp_path, domain_size=6,
        )

"""Keep the published two-sea OSWEC motion gate in ordinary Python CI."""

from hashlib import sha256
from pathlib import Path

import numpy as np

from wecsim import MultiPMWave, PMWave, WEC, WorldPoint
from wecsim.environment.irregularWave import (
    IrregularComponents, synthesize_multiple_irregular_response,
)


FIXTURE = Path(__file__).parent / "fixtures/matlab/oswec_multi_wave"
HYDRO = FIXTURE / "oswec.h5"
SOURCE = FIXTURE / "source.npz"
HYDRO_SHA256 = "b9ccf3a32e0893c23b1645c15a2baef4a0c6f678b90d177d5797e3c65de4b7d3"
SOURCE_SHA256 = "43c9afd588beeb1a837f0f670f33e4bf8d88340b75f2db4c905991a8b7abd7e4"


def _max_error(actual, expected, limit, label):
    error = float(np.max(np.abs(actual - expected)))
    assert error < limit, f"{label}: {error} >= {limit}"


def test_published_two_sea_oswec_motion_and_pto(tmp_path):
    assert sha256(HYDRO.read_bytes()).hexdigest() == HYDRO_SHA256
    assert sha256(SOURCE.read_bytes()).hexdigest() == SOURCE_SHA256
    with np.load(SOURCE) as saved:
        data = {name: saved[name] for name in saved.files}
    np.testing.assert_array_equal(
        data["settings"], [[12700, 1.85e6, 1.85e6, 1.85e6,
                            1, .01, 120000, 0]],
    )
    time = np.arange(1001) * .1
    for name in ("body1", "body2", "pto", "wave1_elevation", "wave2_elevation"):
        np.testing.assert_allclose(data[name][:, 0], time, rtol=0, atol=1e-10)

    seas = []
    components = []
    for number, height, direction in ((1, 2, 0), (2, 1, 90)):
        values = data[f"wave{number}_components"]
        assert values.shape == (500, 4)
        phase_file = tmp_path / f"phase{number}.csv"
        np.savetxt(phase_file, values[:, 3], delimiter=",")
        seas.append(PMWave(height, 3, direction=direction,
                           phase_file=phase_file.name))
        components.append(IrregularComponents(
            omega=values[:, 0], spectral_amplitude=values[:, 1],
            d_omega=values[:, 2], directions=np.array([direction]),
            spreading=np.array([1.]), phase=values[:, 3, None],
        ))

    wec = WEC("Published two-sea OSWEC")
    flap = wec.body("flap", HYDRO, mass=12700, inertia=(1.85e6,) * 3)
    base = wec.body("base", HYDRO, mass=999, inertia=(999,) * 3)
    wec.fixed_hinge(
        flap, base, location=WorldPoint(0, 0, -10),
        pto_location=WorldPoint(0, 0, -8.9), damping=120000,
        added_mass_scheme="simulink_delay",
    )
    result = wec.run(
        MultiPMWave(tuple(seas), excitation_interpolation="spline"),
        dt=.1, end_time=100, ramp_time=.1, radiation_memory=40,
        base_dir=tmp_path,
    )
    np.testing.assert_allclose(result.time, time, rtol=0, atol=1e-10)
    expected_elevation = (data["wave1_elevation"][:, 1]
                          + data["wave2_elevation"][:, 1])
    _max_error(result.wave_elevation, expected_elevation, 1e-11,
               "combined wave elevation (m)")

    # Source body-force columns are independent of the motion integration.
    for number in (1, 2):
        incident = synthesize_multiple_irregular_response(
            HYDRO, tuple(components), dt=.1, end_time=100, ramp_time=.1,
            body_number=number, excitation_interpolation="spline",
        )
        _max_error(incident.excitation_force, data[f"body{number}"][:, 19:25],
                   1e-6, f"body {number} excitation (N or N m)")

    source_flap, source_base, source_pto = (
        data["body1"], data["body2"], data["pto"]
    )
    for dof, position_limit, velocity_limit in (
        (0, .012, .04), (2, .0014, .004), (4, .0025, .008),
    ):
        _max_error(result.bodies["flap"].position[:, dof],
                   source_flap[:, 1 + dof], position_limit,
                   f"flap position {dof}")
        _max_error(result.bodies["flap"].velocity[:, dof],
                   source_flap[:, 7 + dof], velocity_limit,
                   f"flap velocity {dof}")
    for dof in (1, 3, 5):
        np.testing.assert_allclose(result.bodies["flap"].position[:, dof],
                                   source_flap[:, 1 + dof], rtol=0, atol=1e-10)
        np.testing.assert_allclose(result.bodies["flap"].velocity[:, dof],
                                   source_flap[:, 7 + dof], rtol=0, atol=1e-10)
    np.testing.assert_allclose(result.bodies["base"].position,
                               source_base[:, 1:7], rtol=0, atol=1e-10)
    np.testing.assert_allclose(result.bodies["base"].velocity,
                               source_base[:, 7:13], rtol=0, atol=1e-10)
    _max_error(result.ptos["hinge"].force, source_pto[:, 17], 900,
               "PTO torque (N m)")
    _max_error(-result.ptos["hinge"].absorbed_power, source_pto[:, 23], 350,
               "PTO mechanical power (W)")

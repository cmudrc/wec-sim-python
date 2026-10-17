"""Separate ode45 MaxStep from the published nonlinear-force sample clock."""

import json
import os
from pathlib import Path

import numpy as np
import pytest

from wecsim import RegularCICWave, RegularWave, WEC, WorldPoint
from wecsim.hydrodynamics.nonlinearHydro import HeaveMeshHydro


ROOT = Path(__file__).resolve().parents[1]
SOURCE = os.environ.get("WEC_SIM_ELLIPSOID_MAXSTEP_DIR")
APPLICATIONS = os.environ.get("WEC_SIM_APPLICATIONS_DIR")
SELECTED_CASE = os.environ.get("WEC_SIM_ELLIPSOID_CASE")
REFINED = ROOT / "tests/fixtures/matlab/ellipsoid_ode45_maxstep_only/refined_005_common_grid.npz"
pytestmark = pytest.mark.skipif(
    not (SOURCE and APPLICATIONS), reason="fresh pinned MATLAB max-step audit absent",
)


def _read(folder, stem, shape):
    values = np.loadtxt(folder / f"{stem}.csv", delimiter=",", ndmin=2)
    assert values.shape == shape and np.isfinite(values).all(), stem
    return values


@pytest.mark.parametrize(
    "wave_case", [SELECTED_CASE] if SELECTED_CASE else ["Regular", "RegularCIC"],
)
def test_solver_max_step_alone_does_not_close_ellipsoid_gap(wave_case):
    folder = Path(SOURCE) / wave_case
    data = {}
    for mode in ("published", "maxstep005"):
        data[mode] = {
            "body": _read(folder, f"{mode}_body", (3001, 55)),
            "pto": _read(folder, f"{mode}_pto", (3001, 25)),
            "wave": _read(folder, f"{mode}_wave", (3001, 2)),
            "settings": _read(folder, f"{mode}_settings", (1, 5))[0],
        }
        np.testing.assert_allclose(
            data[mode]["settings"],
            [.05, .05, .05, .05, .05 if mode == "published" else .005],
            rtol=0, atol=1e-12,
        )
    published, maxstep = data["published"], data["maxstep005"]
    time = np.arange(3001) * .05
    for mode in (published, maxstep):
        for name in ("body", "pto", "wave"):
            np.testing.assert_allclose(mode[name][:, 0], time, rtol=0, atol=1e-10)
    np.testing.assert_allclose(
        published["wave"][:, 1], maxstep["wave"][:, 1],
        rtol=0, atol=1e-11,
    )

    key = "regular" if wave_case == "Regular" else "regular_cic"
    with np.load(REFINED, allow_pickle=False) as archive:
        fine = {field: archive[f"{key}_{field}"] for field in (
            "time", "heave", "speed", "restoring", "pto_force", "wave")}
    np.testing.assert_allclose(fine["time"], time, rtol=0, atol=1e-10)
    np.testing.assert_allclose(published["wave"][:, 1], fine["wave"],
                               rtol=0, atol=1e-11)

    app = Path(APPLICATIONS) / "Nonlinear_Hydro"
    mesh = HeaveMeshHydro.from_stl(
        app / "geometry/elipsoid.stl", center_z=-2, rho=1025,
        gravity=9.81, depth=70, period=6, height=4, ramp_time=50,
        mass=None, drag_coefficient=1, drag_area=np.pi * 25,
    )
    for mode in (published, maxstep):
        body = mode["body"]
        samples = np.arange(1, len(time), 10)
        preceding = np.array([
            mesh.forces(body[i - 1, 0], body[i - 1, 3] + 2,
                        body[i - 1, 9])[0]
            for i in samples
        ])
        current = np.array([
            mesh.forces(body[i, 0], body[i, 3] + 2, body[i, 9])[0]
            for i in samples
        ])
        assert np.max(np.abs(preceding + body[samples, 39])) < 1e-6
        assert np.max(np.abs(current + body[samples, 39])) > 1_000

    wec = WEC("ellipsoid ode45 max-step isolation")
    body = wec.body(
        "ellipsoid", app / "hydroData/ellipsoid.h5", mass="equilibrium",
        inertia=(1.375264e6, 1.375264e6, 1.341721e6),
        geometry_file=app / "geometry/elipsoid.stl",
        nonlinear_hydro="instantaneous", drag_coefficient=1,
        drag_area=np.pi * 25,
    )
    wec.coordinate("heave", body.move("heave"))
    wec.pto("PTO1", WorldPoint(0, 0, -12.5), body.at(0, 0, 0),
            damping=1_200_000)
    wave = RegularWave(4, 6) if wave_case == "Regular" else RegularCICWave(4, 6)
    python = wec.run(
        wave, dt=.05, end_time=150, ramp_time=50, rho=1025,
        radiation_memory=60 if wave_case == "RegularCIC" else None,
    )
    position = python.bodies["ellipsoid"].position[:, 2]
    speed = python.bodies["ellipsoid"].velocity[:, 2]
    force = python.ptos["PTO1"].force
    metrics = {
        "source_published_to_fine_heave": float(np.max(np.abs(
            published["body"][:, 3] - fine["heave"]))),
        "source_maxstep_to_fine_heave": float(np.max(np.abs(
            maxstep["body"][:, 3] - fine["heave"]))),
        "source_published_to_maxstep_heave": float(np.max(np.abs(
            published["body"][:, 3] - maxstep["body"][:, 3]))),
        "python_to_fine_heave": float(np.max(np.abs(position - fine["heave"]))),
        "python_to_published_heave": float(np.max(np.abs(
            position - published["body"][:, 3]))),
        "python_to_maxstep_heave": float(np.max(np.abs(
            position - maxstep["body"][:, 3]))),
        "python_to_fine_speed": float(np.max(np.abs(speed - fine["speed"]))),
        "python_to_fine_pto_force": float(np.max(np.abs(force - fine["pto_force"]))),
    }
    print(f"ellipsoid {wave_case} max-step isolation: {json.dumps(metrics, sort_keys=True)}")

    # These diagnostic limits were fixed before the RegularCIC MaxStep run.
    # They demonstrate a remaining sampled-force gap, not published parity.
    assert metrics["source_published_to_maxstep_heave"] > .003
    assert metrics["source_published_to_fine_heave"] > .012
    assert metrics["source_maxstep_to_fine_heave"] > .012
    assert metrics["python_to_fine_heave"] < (.003 if wave_case == "Regular" else .0035)
    assert metrics["python_to_fine_speed"] < (.004 if wave_case == "Regular" else .0045)
    assert metrics["python_to_fine_pto_force"] < (5_000 if wave_case == "Regular" else 6_000)
    assert metrics["python_to_published_heave"] > 3 * metrics["python_to_fine_heave"]
    assert metrics["python_to_maxstep_heave"] > 3 * metrics["python_to_fine_heave"]

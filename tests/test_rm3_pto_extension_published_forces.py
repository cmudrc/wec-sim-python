"""Check published RM3 free-decay radiation and added-mass forces."""

from pathlib import Path

import numpy as np
import pytest

from wecsim.hydrodynamics.bodyClass import BodyClass


ROOT = Path(__file__).parent
FORCES = (ROOT / "fixtures/matlab"
          / "rm3_pto_extension_published_forces.npz")
MOTION = ROOT / "fixtures/matlab/rm3_pto_extension_step.npz"
# This checked-in HDF5 has the same active noB2B added mass and IRF as the
# pinned BEMIO RM3 file. The paired force gate also catches a changed input.
HYDRO = ROOT.parent / "examples/data/rm3_end_stops.h5"
DT = .1


def _hydro_force_law(body_number):
    body = BodyClass(str(HYDRO))
    body.bodyNumber = body_number
    body.bodyTotal = 2
    body.readH5file()
    body.mass = "equilibrium"
    body.hydroStiffness = np.zeros((6, 6))
    body.viscDrag = {
        "Drag": np.zeros((6, 6)), "cd": np.zeros(6),
        "characteristicArea": np.zeros(6),
    }
    body.linearDamping = np.zeros((6, 6))
    time = np.arange(301) * DT
    convolution_time = np.arange(601) * DT
    body.hydroForcePre(
        [], [0], len(convolution_time), convolution_time, [],
        DT, 1000, 9.81, "noWaveCIC",
        np.vstack((time, np.zeros_like(time))),
        body_number, 2, 0, 0, 0,
    )
    kernel = np.asarray(body.hydroForce["irkb"])[:, 2, 2]
    added = np.asarray(body.hydroForce["fAddedMass"])
    assert kernel.shape == (601,) and added.shape == (6, 6)
    return kernel, added


@pytest.mark.parametrize("case,active_body", (("float", 1), ("spar", 2)))
def test_published_free_decay_force_laws(case, active_body):
    with (np.load(FORCES, allow_pickle=False) as saved,
          np.load(MOTION, allow_pickle=False) as archived):
        records = [saved[f"{case}_body{number}"].copy()
                   for number in (1, 2)]
        reference = archived[f"{case}_0p1"]
    assert reference.shape == (301, 9)
    for number, record in enumerate(records, start=1):
        assert record.shape == (301, 6) and np.isfinite(record).all()
        np.testing.assert_allclose(record[:, 0], np.arange(301) * DT,
                                   rtol=0, atol=1e-12)
        # The new force export is from the exact archived published motion.
        np.testing.assert_allclose(record[:, 1:3],
                                   reference[:, 2 * number - 1:2 * number + 1],
                                   rtol=0, atol=1e-10)

    record = records[active_body - 1]
    velocity = record[:, 2]
    acceleration = record[:, 3]
    kernel, added = _hydro_force_law(active_body)

    # The source convolution block uses a trapezoidal IRF history at accepted
    # major steps. The body is constrained to heave in these two no-wave cases.
    radiation = np.empty(len(record))
    for step in range(len(record)):
        radiation[step] = DT / 2 * kernel[0] * velocity[step]
        if step:
            radiation[step] += DT * np.dot(
                kernel[1:step + 1], velocity[step - 1::-1],
            )
    radiation_gap = np.max(np.abs(radiation - record[:, 4]))
    assert radiation_gap < .005, (case, radiation_gap)

    # WEC-Sim shifts twice the translational added-mass trace into rigid
    # inertia; its remaining heave coefficient multiplies a 100 ns delayed
    # acceleration, extrapolated from prior accepted major steps.
    applied = added[2, 2] - 2 * np.trace(added[:3, :3])
    delayed = np.zeros_like(acceleration)
    delayed[1] = acceleration[0]
    factor = 1 - 1e-7 / DT
    delayed[2:] = (acceleration[1:-1]
                   + factor * (acceleration[1:-1] - acceleration[:-2]))
    added_gap = np.max(np.abs(applied * delayed - record[:, 5]))
    assert added_gap < .025, (case, added_gap)

    inactive = records[2 - active_body]
    assert np.max(np.abs(inactive[:, 2])) < 1e-9
    assert np.max(np.abs(inactive[:, 4:6])) < .005

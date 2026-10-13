"""Exercise the public RM3 mooring coupling without a native installation."""

from pathlib import Path

import numpy as np

from wecsim import MoorDyn, NoWave, WEC


HYDRO = Path(__file__).resolve().parents[1] / "examples/data/rm3.h5"


class _ConstantSurgeMooring:
    def __init__(self):
        self.calls = []
        self.closed = False

    def start(self, position, velocity):
        self.initial_position = np.array(position, copy=True)
        self.initial_velocity = np.array(velocity, copy=True)
        return self

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.closed = True

    def step(self, position, velocity, time, dt):
        self.calls.append((np.array(position, copy=True),
                           np.array(velocity, copy=True), time, dt))
        return np.array([100_000., 0, 0, 0, 0, 0])


def test_public_moordyn_load_uses_spar_local_lever_arm(tmp_path):
    library = tmp_path / "moordyn.so"
    lines = tmp_path / "lines.txt"
    library.touch()
    lines.touch()
    final_pitch = {}
    for point_height in (0, 21.5):
        wec = WEC("RM3 mooring lever arm")
        float_body = wec.body("float", HYDRO,
                              inertia=(0, 21_306_090.66, 0))
        spar = wec.body("spar", HYDRO,
                        inertia=(0, 94_407_091.24, 0))
        mooring = MoorDyn(library, lines)
        session = _ConstantSurgeMooring()
        mooring.start = session.start
        mooring.step = session.step
        wec.floating_joint(
            float_body, spar, damping=1_200_000,
            moordyn=mooring, moordyn_point=spar.at(0, 0, point_height),
            radiation_method="convolution",
        )
        result = wec.run(
            NoWave(), dt=.05, end_time=.2, radiation_memory=.2,
            initial_coordinate={"spar_heave": -.21},
        )
        output = dict(result.raw.extra_outputs)
        position = output["moordyn_connection_position"]
        velocity = output["moordyn_connection_velocity"]
        force = output["moordyn_connection_force"]
        assert session.closed
        assert len(session.calls) == 4
        np.testing.assert_allclose(session.initial_position, position[0])
        np.testing.assert_allclose(session.initial_velocity, velocity[0])
        np.testing.assert_allclose(position[0, 2], point_height - 21.5,
                                   rtol=0, atol=1e-12)
        np.testing.assert_array_equal(force[0], np.zeros(6))
        np.testing.assert_array_equal(force[1:, 0], np.full(4, 100_000.))
        for index, (pose, speed, time, dt) in enumerate(session.calls, 1):
            np.testing.assert_allclose(pose, position[index])
            np.testing.assert_allclose(speed, velocity[index])
            np.testing.assert_allclose((time, dt), (result.time[index - 1], .05))
        final_pitch[point_height] = result.bodies["spar"].position[-1, 4]

    # The same horizontal force reverses its pitch moment when the
    # attachment crosses the spar's joint reference height.
    assert final_pitch[0] < -1e-5
    assert final_pitch[21.5] > 1e-5

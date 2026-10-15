"""Keep the local R2025b OWC source-to-source audit distinguishable from parity."""

from hashlib import sha256
from pathlib import Path

import numpy as np


FIXTURES = Path(__file__).parent / "fixtures/matlab"
HYDRO = FIXTURES / "owc_orifice.h5"
LOCAL = FIXTURES / "owc_orifice_local_prefix/source.npz"
ARCHIVED = FIXTURES / "owc_orifice_prefix.npz"


def test_local_source_prefix_has_identical_inputs_but_fails_motion_gates():
    assert sha256(HYDRO.read_bytes()).hexdigest() == (
        "fe702ecf99f10d4cd0bba5d858b588dd37709cd45a5c3727356849e8dc40d676"
    )
    with np.load(LOCAL, allow_pickle=False) as local, \
            np.load(ARCHIVED, allow_pickle=False) as archived:
        assert local["body"].shape == (1201, 13)
        assert local["flex"].shape == (1201, 10)
        assert local["components"].shape == (500, 4)
        np.testing.assert_array_equal(local["components"], archived["components"])
        np.testing.assert_array_equal(local["wave"], archived["wave"])
        np.testing.assert_allclose(local["body"][:, 0], np.arange(1201) * .005,
                                   rtol=0, atol=1e-12)
        np.testing.assert_allclose(
            local["body"][:401, 1:13], archived["body"][:401, 1:13],
            rtol=0, atol=1e-12,
        )

        # These are the limits declared before the local run. The source
        # histories separate despite identical waves, so this is not a
        # passing independent Python/published-MATLAB motion comparison.
        body_position = np.max(np.abs(
            local["body"][:, 1:7] - archived["body"][:, 1:7],
        ))
        body_speed = np.max(np.abs(
            local["body"][:, 7:13] - archived["body"][:, 7:13],
        ))
        flex_position = np.max(np.abs(
            local["flex"][:, 0] - archived["flex"][:, 0],
        ))
        flex_speed = np.max(np.abs(
            local["flex"][:, 1] - archived["flex"][:, 1],
        ))
        assert body_position > 1e-5
        assert body_speed > 1e-5
        assert flex_position > 1e-5
        assert flex_speed > 1e-5

        assert local["archived_orifice"].shape == (1639, 6)
        assert local["orifice"].shape == (1437, 6)
        assert abs(local["archived_orifice"][-1, 1]
                   - local["orifice"][-1, 1]) > .5

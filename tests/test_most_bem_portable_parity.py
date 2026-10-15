"""Run the pinned MATLAB MOST blade-load gate in ordinary Python CI."""

from pathlib import Path

import numpy as np

from wecsim import MostBEM


SOURCE = Path(__file__).parent / "fixtures/matlab/most_bem/source.npz"
MODEL_FIELDS = (
    "radius", "radius_interval", "twist", "chord", "curve_angle",
    "sweep", "curve", "airfoil_index", "airfoil", "hub_radius",
    "tip_radius", "rho_air", "root_tol", "max_iterations",
    "root_func_tol", "root_tol_fallback", "max_fallback_iterations",
    "epsilon",
)


def test_published_most_blade_loads_match_five_matlab_states():
    """Advance all three blades using only archived inputs, not source loads."""
    with np.load(SOURCE, allow_pickle=False) as source:
        model = MostBEM(**{name: source[name].item() if source[name].ndim == 0
                           else source[name] for name in MODEL_FIELDS})
        states = source["hub_state"]
        speeds = source["wind_speeds"]
        pitches = source["blade_pitch"]
        expected = source["root_loads"]
        assert states.shape == (5, 14)
        assert expected.shape == (5, 6, 3)

        for index, (state, speed, pitch) in enumerate(zip(
            states, speeds, pitches, strict=True,
        )):
            if index == 4:
                # The fifth MATLAB state uses the published spatial wind field.
                wind = lambda point: np.array([
                    8 + .002*point[0] + .004*point[1] + .01*(point[2]-150),
                    .1 + .001*point[1], .02*point[2]/150,
                ])
            else:
                wind = np.array([speed, 0, 0])
            actual = model.loads(state, pitch, wind)
            np.testing.assert_allclose(actual, expected[index],
                                       rtol=0, atol=1e-4)

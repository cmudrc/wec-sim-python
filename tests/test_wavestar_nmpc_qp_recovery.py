"""Recovery from the SciPy line-search failure seen in the NMPC source gate."""

import numpy as np
from scipy.optimize import OptimizeResult

import wecsim.wavestarNmpc as nmpc


def test_nmpc_qp_recovers_after_primary_line_search_failure(monkeypatch):
    state = np.array([.05, .4, 0., 0., .3])
    forecast = 3 * np.sin(np.linspace(0, 2 * np.pi, 40))
    expected = nmpc.WaveStarNmpcController().step(state, forecast, 15.)

    original = nmpc.minimize
    methods = []

    def line_search_failure(*args, **kwargs):
        methods.append(kwargs["method"])
        result = original(*args, **kwargs)
        if kwargs["method"] == "L-BFGS-B":
            return OptimizeResult(
                success=False, message="ABNORMAL:",
                x=.5 * result.x,
            )
        return result

    monkeypatch.setattr(nmpc, "minimize", line_search_failure)
    recovered = nmpc.WaveStarNmpcController().step(state, forecast, 15.)
    assert methods == ["L-BFGS-B", "SLSQP"]
    np.testing.assert_allclose(recovered, expected, rtol=0, atol=1e-6)

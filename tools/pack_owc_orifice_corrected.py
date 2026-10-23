"""Pack a corrected OWC MATLAB CI artifact into the checked source fixture."""

import argparse
from pathlib import Path

import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifact_dir", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    names = ("body", "flex", "orifice", "wave", "components")
    source = {
        name: np.loadtxt(args.artifact_dir / f"{name}.csv", delimiter=",")
        for name in names
    }
    source["body"] = source["body"][:, :13]
    time = np.arange(26_001) * .005
    internal_orifice = source["orifice"]
    if internal_orifice.shape[0] > len(time):
        assert internal_orifice.shape[1] == 6
        assert np.all(np.diff(internal_orifice[:, 0]) > 0)
        with np.load(
            Path(__file__).parents[1]
            / "tests/fixtures/matlab/owc_orifice_130s/source.npz",
            allow_pickle=False,
        ) as published:
            orifice_area = float(published["orifice_parameters"][1])
        source["summary"] = np.array([
            np.trapezoid(internal_orifice[:, 5] * 1000,
                         internal_orifice[:, 0]),
            np.max(np.abs(internal_orifice[:, 4] / orifice_area) / 343.2),
        ])
        source["orifice"] = np.column_stack((
            time,
            *[np.interp(time, internal_orifice[:, 0], internal_orifice[:, i])
              for i in range(1, 6)],
        ))
    else:
        source["summary"] = np.loadtxt(
            args.artifact_dir / "summary.csv", delimiter=",",
        )
    assert source["body"].shape == (26_001, 13)
    assert source["flex"].shape == (26_001, 11)
    assert source["orifice"].shape == (26_001, 6)
    assert source["wave"].shape == (26_001, 2)
    assert source["components"].shape == (500, 4)
    assert source["summary"].shape == (2,)
    for name in ("body", "flex", "orifice", "wave"):
        np.testing.assert_allclose(source[name][:, 0], time, rtol=0, atol=1e-10)
        assert np.isfinite(source[name]).all()
    assert np.isfinite(source["components"]).all()
    assert source["summary"][0] > 0
    assert 0 <= source["summary"][1] < .3
    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.output, **source)


if __name__ == "__main__":
    main()

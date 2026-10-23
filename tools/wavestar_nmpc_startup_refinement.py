"""Probe independent physical step refinement before WaveStar NMPC feedback.

The checked source supplies only its hydrodynamic file and wave phases to the
Python runs. Source motion and requests are loaded after each run for comparison.
"""

import argparse
import json
import platform
import sys
from hashlib import sha256
from pathlib import Path

import numpy as np
import scipy

from wecsim.dynamics.wavestar import run_wavestar_published
from wecsim.environment.irregularWave import jonswap_equal_energy_components


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tests"))
from test_wavestar_split_step_portable import _DelayedSampledPTO  # noqa: E402


FIXTURE = (ROOT / "tests/fixtures/matlab/"
           "wavestar_nmpc_reference_convergence_225s/"
           "source_0p25ms_compatible_amd7763.npz")
HYDRO = (ROOT / "tests/fixtures/matlab/"
         "wavestar_nmpc_split_step_225s/wavestar.h5")
SOURCE_SHA256 = "48c02a9396fd166d809abd8ce4cdd136934f9d819c6a7b9c73f01703907353b6"
HYDRO_SHA256 = "b2541ff01c319987d9157d1be9fbc9ec39a7660d0572d8e316c428b8aa3359d1"
STEPS = (.001, .0005, .00025, .000125, .0000625)
END_TIME = 15.


def audit():
    if sha256(FIXTURE.read_bytes()).hexdigest() != SOURCE_SHA256:
        raise ValueError("The checked WaveStar source fixture changed")
    if sha256(HYDRO.read_bytes()).hexdigest() != HYDRO_SHA256:
        raise ValueError("The WaveStar hydrodynamic fixture changed")
    with np.load(FIXTURE, allow_pickle=False) as fixture:
        source = {key: fixture[key].copy() for key in fixture.files}
    sea = jonswap_equal_energy_components(
        HYDRO, significant_height=.1042, peak_period=1.836,
        directions=np.array([0.]), spreading=np.array([1.]),
        gamma=3.3, phase=source["components"][:, 3, None],
    )
    count = round(END_TIME / .05) + 1
    np.testing.assert_allclose(source["time"][:count],
                               np.arange(count) * .05, rtol=0, atol=1e-9)
    results = {}
    for step in STEPS:
        controller = _DelayedSampledPTO(plant_dt=step)
        response = run_wavestar_published(
            HYDRO, sea, dt=step, end_time=END_TIME, ramp_time=25.,
            g=9.80665, pto_controller=controller,
            output_stride=round(.05 / step),
        )
        np.testing.assert_allclose(response.time, source["time"][:count],
                                   rtol=0, atol=1e-9)
        actual = {
            "pitch": response.angle,
            "pitch_speed": response.angular_speed,
            "stroke": response.pto_stroke,
            "axial_force": response.pto_force,
            "command": np.asarray(controller.commands),
        }
        results[f"{step:g}"] = {
            name: {
                "gap_at_15s": float((values - source[name][:count])[-1]),
                "max_abs_gap_to_15s": float(np.max(np.abs(
                    values - source[name][:count],
                ))),
            }
            for name, values in actual.items()
        }
    return {
        "source_sha256": SOURCE_SHA256,
        "hydro_sha256": HYDRO_SHA256,
        "source_plant_step_s": .00025,
        "python_plant_steps_s": STEPS,
        "startup_end_s": END_TIME,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "results": results,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = audit()
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()

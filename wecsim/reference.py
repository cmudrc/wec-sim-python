"""Command-line runner for the three supported reference model families.

Run ``python -m wecsim.reference MODEL --h5 FILE --output FILE``.
The CSV and adjacent JSON record numeric results and reproducibility inputs.
These focused models are distinct from WEC-Sim's general Simulink runner.
"""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess

import numpy as np

from .hingePitch import solve_hinged_pitch_from_excitation
from .irregularWave import pm_equal_energy_components, synthesize_irregular_response
from .linearHeave import solve_heave_free_decay
from .rm3Regular import solve_rm3_regular


def _sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _git_state():
    root = Path(__file__).resolve().parents[1]
    try:
        revision = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=root,
            capture_output=True, text=True, check=False,
        )
        changes = subprocess.run(
            ["git", "status", "--porcelain"], cwd=root,
            capture_output=True, text=True, check=False,
        )
    except OSError:
        return None, None
    commit = revision.stdout.strip() if revision.returncode == 0 else None
    dirty = bool(changes.stdout) if changes.returncode == 0 else None
    return commit, dirty


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="model", required=True)
    for name in ("rm3", "oswec", "sphere"):
        sub = subparsers.add_parser(name)
        sub.add_argument("--h5", type=Path, required=True, help="hydrodynamic HDF5 file")
        sub.add_argument("--output", type=Path, required=True, help="result CSV path")
    subparsers.choices["rm3"].add_argument("--end-time", type=float, default=400.0)
    subparsers.choices["rm3"].add_argument(
        "--b2b", action="store_true", help="enable RM3 body-to-body hydrodynamic coupling",
    )
    subparsers.choices["oswec"].add_argument("--end-time", type=float, default=400.0)
    subparsers.choices["oswec"].add_argument("--seed", type=int, default=7)
    subparsers.choices["sphere"].add_argument("--end-time", type=float, default=40.0)
    subparsers.choices["sphere"].add_argument(
        "--initial-displacement", type=float, required=True,
    )
    args = parser.parse_args(argv)
    h5_file = args.h5.expanduser().resolve(strict=True)
    output = args.output.expanduser().resolve()
    if output == h5_file:
        parser.error("CSV output cannot overwrite the hydrodynamic input")
    metadata_path = output.with_suffix(".json")
    if metadata_path == output or metadata_path.resolve() == h5_file:
        metadata_path = output.with_suffix(".provenance.json")
    if metadata_path == output or metadata_path.resolve() == h5_file:
        parser.error("provenance output cannot overwrite the hydrodynamic input")
    output.parent.mkdir(parents=True, exist_ok=True)

    if args.model == "rm3":
        response = solve_rm3_regular(h5_file, end_time=args.end_time, b2b=args.b2b)
        columns = ["time"]
        arrays = [response.time]
        for name, values in (("position", response.body_position),
                             ("velocity", response.body_velocity)):
            for body in range(2):
                for dof, label in ((0, "surge"), (2, "heave"), (4, "pitch")):
                    columns.append(f"body{body + 1}_{label}_{name}")
                    arrays.append(values[:, body, dof])
        columns.append("pto_heave_force_on_body1")
        arrays.append(response.pto_force)
        parameters = {"end_time": args.end_time, "dt": 0.1,
                      "wave_height": 2.5, "wave_period": 8.0,
                      "ramp_time": 100.0, "pto_damping": 1_200_000.0,
                      "b2b": args.b2b}
    elif args.model == "oswec":
        components = pm_equal_energy_components(
            h5_file, significant_height=2.5, peak_period=8.0,
            directions=[0, 30, 90], spreading=[0.1, 0.2, 0.7], seed=args.seed,
        )
        wave = synthesize_irregular_response(
            h5_file, components, dt=0.1, end_time=args.end_time, ramp_time=100.0,
        )
        response = solve_hinged_pitch_from_excitation(
            h5_file, wave.excitation_force, hinge_z=-8.9,
            body_mass=127_000, pitch_inertia=1.85e6, pto_damping=12_000,
        )
        columns = ["time", "wave_elevation", "surge_position", "heave_position",
                   "pitch_angle", "pitch_velocity", "pto_pitch_torque"]
        arrays = [response.time, wave.elevation, response.center_position[:, 0],
                  response.center_position[:, 2], response.angle,
                  response.angular_velocity, response.pto_torque]
        parameters = {"end_time": args.end_time, "dt": 0.1,
                      "significant_height": 2.5, "peak_period": 8.0,
                      "directions": [0, 30, 90], "spreading": [0.1, 0.2, 0.7],
                      "ramp_time": 100.0, "seed": args.seed,
                      "pto_damping": 12_000}
    else:
        response = solve_heave_free_decay(
            h5_file, args.initial_displacement, end_time=args.end_time,
        )
        columns = ["time", "heave_position", "heave_velocity", "total_force"]
        arrays = [response.time, response.position, response.velocity,
                  response.force_total]
        parameters = {"end_time": args.end_time, "dt": 0.01,
                      "initial_displacement": args.initial_displacement,
                      "cic_end_time": 15.0}

    values = np.column_stack(arrays)
    if not np.isfinite(values).all():
        raise RuntimeError("the reference case produced nonfinite output")
    np.savetxt(output, values, delimiter=",", header=",".join(columns),
               comments="", fmt="%.17g")
    commit, dirty = _git_state()
    metadata = {
        "model": args.model,
        "h5_file": str(h5_file),
        "h5_sha256": _sha256(h5_file),
        "git_commit": commit,
        "git_dirty": dirty,
        "numpy_version": np.__version__,
        "parameters": parameters,
        "csv_columns": columns,
        "time_steps": len(values),
    }
    metadata_path.write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8",
    )
    print(output)


if __name__ == "__main__":
    main()

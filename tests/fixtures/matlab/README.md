# RM3 regular-wave MATLAB trajectory

`rm3_regular_wave.npz` contains selected motion and PTO channels at every
sample of the 400 s published RM3 regular-wave baseline. It was extracted
without rounding from the successful
[pinned R2025b MATLAB run](https://github.com/cmudrc/wec-sim-python/actions/runs/38212367776)
(`matlab-RM3-baseline` artifact, run commit `dbcd62bfeba7bcc857f699522d17c06add834c39`).
The workflow checks out WEC-Sim core
`0753b2e47f2457c078751dcfe5d251d1767b80ab`. Inputs are a 2.5 m,
8 s regular wave, 0.1 s output step, 100 s ramp, and a 1.2 MN s/m linear PTO.

The archive contains `time`, `body_position`, `body_velocity`,
`body_heave_excitation_force`, `pto_stroke`, `pto_speed`, `pto_force`, and
`pto_power`. Body arrays are ordered float/spar, then time; position and
velocity add surge/sway/heave/roll/pitch/yaw as the last axis. These are
source CSV columns 0, 1–6, 7–12, and 21 for bodies, and columns 3, 9, 15,
and 21 for the PTO.
The source files in the artifact have these SHA-256 hashes:

| Source file | SHA-256 |
| --- | --- |
| `RM3_RM3_body1.csv` | `f3f31922d9fac06b13b1ad6393443256084d30c86ffe1788eca4a97e3273ff50` |
| `RM3_RM3_body2.csv` | `60594f5c7cebf17e3dd8347210abab4b1324b3e3ba50b9cc8631c7b66d88ecaf` |
| `RM3_RM3_pto1.csv` | `d67f8833edb5cd67bdba5a48b76a6a440cd7eed68dce4800a108d478a014f600` |

To regenerate, load each CSV with `numpy.loadtxt(..., delimiter=",")`, check
that all three time columns match, and save those columns using
`numpy.savez_compressed` under the names and ordering above. The archive SHA-256
is `bc8db627eb710b002e6ed91b0a352999ea62a56c09baf5e859677fa1f1442ba5`.

The always-on test uses the bundled `examples/data/rm3.h5` (SHA-256
`3bdfdcfa3dcd2db1851ab2026eb6b4fdc565ad86be8503e09b42231feee429fb`).
That HDF5 is not byte-identical to the current WEC-Sim example; its regular-wave
hydrodynamic data still reproduce this source trajectory within the stated
per-channel gates. The live paired tests continue to use pinned MATLAB inputs.

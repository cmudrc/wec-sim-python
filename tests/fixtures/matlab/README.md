# Portable MATLAB reference data

## OSWEC nonlinear visualization motion

`oswec_nonlinear_viz_published_motion.npz` keeps all 1,201 samples of the
published 120 s, 0.1 s `ode4` flap motion and its principal force channels.
It was copied without rounding from the untouched R2025b
[published-case run](https://github.com/cmudrc/wec-sim-python/actions/runs/38259283403)
(`matlab-oswec-nonlinear-viz-pressure` artifact, run commit
`1083a1f29725785cbba0cb67f8e389a9cade16fc`). The workflow pins WEC-Sim
core `0753b2e47f2457c078751dcfe5d251d1767b80ab` and Applications
`d53d4d4c9eda2581f04204f5d394a6ef84bb099e`. The source `source.mat`
hash is `f5c2944777f35741a48034e64fb4d3d7ec165466f6f2fc76f7be1f38383ad466`.
The archive contains `time`, `pose`, `velocity`, `acceleration`,
`forceExcitation`, `forceRadiationDamping`, `forceAddedMass`, and
`forceRestoring` for the six flap degrees of freedom. Its SHA-256 is
`63f3669b192ab271a92f5ea9c4650d1d183ef68f0abc4a6b0faedcd2b4948e4f`.
The refined-step source job gates its 0.1 s run against these published
channels before evaluating step sensitivity or independent Python motion.

## RM3 regular wave

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

## Sphere free decay

`sphere_free_decay.npz` retains every sample of the five published Sphere
`noWaveCIC` free-decay cases, ordered `0m`, `1m`, `1m-ME`, `3m`, `5m`.
The successful [pinned R2025b run](https://github.com/cmudrc/wec-sim-python/actions/runs/38196174533)
(`matlab-Sphere-baseline` artifact, run commit
`33422a97e5e1273c2828b923d2f96d7c37f31caf`) used WEC-Sim core
`0753b2e47f2457c078751dcfe5d251d1767b80ab` and Applications
`d53d4d4c9eda2581f04204f5d394a6ef84bb099e`. Its BEMIO-generated
`sphere.h5` is copied to `examples/data/sphere.h5` without modification
(SHA-256 `7287c5b02c74616e385f0d582df6a4679d7a4ceb93d8c663833a944838598e06`).

The archive stores `case_names`, a common `time`, and arrays ordered by case,
time, then WEC-Sim degree of freedom for `body_position` and `body_velocity`.
`body_total_heave_force` has case and time axes. Body CSV columns 0, 1–6,
7–12, and 15 were copied as float64 without rounding. Source hashes are:

| Source file | SHA-256 |
| --- | --- |
| `Sphere_0m_body1.csv` | `dd2f2f5d793c769df4ffe9610d4ad5773bf7d4780dab1b224babe88e10ba24ed` |
| `Sphere_1m_body1.csv` | `8eeacfa184077705ec88b8a87f35c32eff097c2a0ee344bca1216a0406ebf427` |
| `Sphere_1m-ME_body1.csv` | `8eeacfa184077705ec88b8a87f35c32eff097c2a0ee344bca1216a0406ebf427` |
| `Sphere_3m_body1.csv` | `4205f502b2c4c0fd14d14c1ba2c01a35bab5331ede5f1d36cb8cad9eb968e44b` |
| `Sphere_5m_body1.csv` | `6d66112742795e9d57009e08b1e88c5f68040669d5fdd589f7d962840d438876` |

The archive SHA-256 is
`95057bc8d58b1e36243eefb688a184774a557dc63542128efff3eac732095154`.
The two 1 m source trajectories are identical: the published Morison element
has surge-only coefficients and does not change heave free decay.

## OSWEC directional irregular wave

`oswec_irregular_wave.npz` keeps all 4,001 samples of the published 400 s
OSWEC PM sea, both body motions, the flap's six excitation forces, and the
hinge PTO. The source is the successful [pinned R2025b run](https://github.com/cmudrc/wec-sim-python/actions/runs/37932186068)
(`matlab-WAVE_MARKERS_OSWEC-baseline` artifact, run commit
`9872f98ccd65aa57e49f0d7dbd5db8919d87773c`). It checked out WEC-Sim
core `0753b2e47f2457c078751dcfe5d251d1767b80ab`. The published input
has three incident directions and no phase seed, so the archive retains its
500-by-3 phase realization. `examples/data/oswec.h5` is an unchanged copy of
that pinned core revision's `examples/OSWEC/hydroData/oswec.h5` (SHA-256
`28aa028d404a88f877924bc435e7ad3a4c50e26f9706c1f39ef7324e5ec60750`).

The archive stores `time`, `wave_elevation`, `wave_components` (frequency,
amplitude, width, then three phases), `wave_directions` (direction and
spreading), `body_position` and `body_velocity` (flap/base, time, six degrees
of freedom), `flap_excitation_force`, and PTO stroke, speed, force, and power.
The body CSV uses columns 0, 1–6, 7–12, and 19–24; the PTO uses columns 0,
5, 11, 17, and 23. All source CSV values were loaded as float64 without
rounding. Their SHA-256 hashes are:

| Source file | SHA-256 |
| --- | --- |
| `OSWEC_OSWEC_body1.csv` | `9b005121c637987aff76915f1b0295200970c7e3bec3711956e9cc25b8570edf` |
| `OSWEC_OSWEC_body2.csv` | `cc8d525d09cc876ae50766f6d1067face24f9a10a823afa19a06d702e99eb33b` |
| `OSWEC_OSWEC_pto1.csv` | `3a28e1d7dbcff757fb980f35ce828f895fb0b9f1481ca0f22ceb90e149775766` |
| `OSWEC_wave_components.csv` | `29a5ce24a901c819953da8ed365bc4c22e8b2beb66f2a8529124db6c19191b00` |
| `OSWEC_wave_directions.csv` | `fb581ce2fe25dcc31394e43c243ff2281627613fe0df5a543f712d4e9343cebf` |
| `OSWEC_wave_elevation.csv` | `35aae2f59b1b66aa2ab6107bca0eca712f2ffbb8c58bce208d2f9d9fcbfdc62b` |

The archive SHA-256 is
`070365a8fd15dddab8b4fe89a6000eb8572f291ee6f2068d3520e13612424b01`.

# WaveStar auxiliary-mass diagnostic

`source.npz` is a local MATLAB R2025b Update 7 run of the pinned
`WECCCOMP_Nonlinear_Model_Predictive` application through 9.95 s. The
[probe](../../../matlab/wavestar_auxiliary_mass_probe.m) copies the application
to a fresh directory and keeps its 50 ms `ode8` solver, sea-state settings,
phase seed, hydrodynamic HDF5, float, arm, constraints, and controller.
It changes only the BC rod and motor masses and three inertia entries from
`0.0001` to `1e-9` in source units. The plotting script is disabled and
geometry filename case is fixed for the local filesystem.

The same probe accepts `20` as a fifth argument for a bounded active-control
experiment. `active_20s.npz` stores its 401 pitch, pitch-speed, and torque
request samples along with the unchanged sea components and wave history.
Generate it from the probe's `aux_mass_probe.mat` with `scipy.io.loadmat`,
selecting `time`, `pitch`, `speed`, `command`, `components`, and `wave` and
flattening the first four vectors. The ordinary 9.95 s fixture remains the
default probe duration. The 20 s archive SHA-256 is
`bd5ef7e1972a0b3c0e009194d2bb7dbb781e87a3dad19f28e0d32dfe2d4396ba`.

Pinned revisions: WEC-Sim core `0753b2e47f2457c078751dcfe5d251d1767b80ab`,
WEC-Sim_Applications `d53d4d4c9eda2581f04204f5d394a6ef84bb099e`.
The HDF5 is `../wavestar_fixed_input_step/wavestar.h5`. Run the probe with
`coreDir`, `applicationsDir`, `repoDir`, and a fresh `outDir`, then convert
`aux_mass_probe.mat` to the compressed NumPy archive:

```python
from scipy.io import loadmat
import numpy as np

data = loadmat("<outDir>/aux_mass_probe.mat")
np.savez_compressed(
    "source.npz", time=data["time"].ravel(), pitch=data["pitch"].ravel(),
    speed=data["speed"].ravel(), acceleration=data["acceleration"],
    auxiliary_mass=np.float64(1e-9),
)
```

| 9.95 s diagnostic | Published source | Near-zero auxiliaries | Python source-ode8 |
| --- | ---: | ---: | ---: |
| Initial float pitch acceleration (rad/s²) | 0.3383701634 | 0.3384090166 | 0.3384090170 |
| Maximum pitch difference from Python (rad) | 1.42056e-5 | 1.44994e-5 | — |
| Maximum pitch-speed difference from Python (rad/s) | 1.18577e-4 | 1.18213e-4 | — |

Changing the auxiliaries shifts MATLAB's initial acceleration by
`3.88532e-5` rad/s², almost exactly the original initial Python/source
gap. The probe changes gravity and inertia together: at the initial linkage
pose, the projected auxiliary mass is `3.3634463e-6` kg m² and the gravity
load is `1.28195481e-4` N m. Dividing by the base generalized mass
`3.32874780505` kg m² attributes about `3.42e-7` rad/s² of the acceleration
change to inertia and `3.85e-5` rad/s² to gravity. Combining both predicts
`0.33837016345` rad/s², within `5.5e-11` of the published source.

Across the full prefix, changing both auxiliaries shifts MATLAB pitch by at
most `2.48887e-6` rad and pitch speed by at most `1.05876e-5` rad/s. The
Python/source gap does not close with the old piecewise-linear solver-stage
wave forcing. This diagnostic isolates the omitted bodies as the cause of
the initial derivative gap; it does not establish full 225 s NMPC parity.

In the 20 s variant, the modified MATLAB run and the unchanged published run
have identical components and wave samples. Their pitch differs by at most
`2.48887e-6` rad through 15 s, yet their controller requests differ by
`0.03235` N m at 16.55 s and up to `0.17171` N m through 20 s. Their pitch
separates by up to `0.00116142` rad by 20 s. This source-only perturbation
shows that small physical changes can strongly alter the early NMPC request
history. On the unchanged published source, the corrected independent Python
`wecsim_ode8` path has a maximum 20 s pitch gap of `0.00064794` rad, versus
`0.00118738` rad before the auxiliary-body and stage-force corrections. This
bounded result is not a gate for the 225 s published input, and it does not
justify changing that gate or the case badge.

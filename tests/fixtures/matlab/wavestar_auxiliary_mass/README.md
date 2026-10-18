# WaveStar auxiliary-mass diagnostic

`source.npz` is a local MATLAB R2025b Update 7 run of the pinned
`WECCCOMP_Nonlinear_Model_Predictive` application through 9.95 s. The
[probe](../../../matlab/wavestar_auxiliary_mass_probe.m) copies the application
to a fresh directory and keeps its 50 ms `ode8` solver, sea-state settings,
phase seed, hydrodynamic HDF5, float, arm, constraints, and controller.
It changes only the BC rod and motor masses and three inertia entries from
`0.0001` to `1e-9` in source units. The plotting script is disabled and
geometry filename case is fixed for the local filesystem.

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
gap. Across the full prefix it changes MATLAB pitch by at most `2.48887e-6`
rad and pitch speed by at most `1.05876e-5` rad/s. The Python/source pitch
gap does not close. This rules out omitted BC/motor inertia as the leading
9.95 s trajectory discrepancy; it does not establish full 225 s NMPC parity.

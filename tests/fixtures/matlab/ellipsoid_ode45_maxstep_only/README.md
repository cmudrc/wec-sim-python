# Ellipsoid ode45 MaxStep isolation

`refined_005_common_grid.npz` retains the 0.05 s common-grid time, heave,
heave speed, restoring force, PTO force, and wave elevation from the full
150 s, `simu.dt=0.005` R2025b [Regular and RegularCIC source
artifacts](https://github.com/cmudrc/wec-sim-python/actions/runs/38110743157).
Each array is every tenth row of the corresponding `step005_*.csv` export.
The fixture SHA-256 is
`7a80b56e69042e230108907fb8a90af062380cfedd1c29fcaf505e029d4446db`.

The fresh workflow checks out WEC-Sim core
`0753b2e47f2457c078751dcfe5d251d1767b80ab` and Applications
`d53d4d4c9eda2581f04204f5d394a6ef84bb099e`. It generates one shared
`ellipsoid.h5` from the pinned WAMIT output, then runs the original 0.05 s
`ode45` input twice. The second run changes only Simulink `MaxStep` from
0.05 to 0.005 s after WEC-Sim initialization; output, nonlinear-force, and
radiation-convolution intervals remain 0.05 s. Saved settings and all 3,001
source samples are checked by `tests/test_ellipsoid_ode45_maxstep_only.py`.

The test also advances the public Python model from the pinned BEM and wave
input. Its limits are a diagnostic of solver and force-sample timing. A
passing result does **not** count either published 0.05 s input as a parity
gate, because the source's restoring force uses the preceding sample.

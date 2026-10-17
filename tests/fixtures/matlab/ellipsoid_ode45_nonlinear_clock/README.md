# Ellipsoid source-clock reference

`source_maxstep_common_grid.npz` contains the 3,001 common 0.05 s samples
from the pinned R2025b `Regular` and `RegularCIC` jobs in [run
39054809363](https://github.com/cmudrc/wec-sim-python/actions/runs/39054809363).
Each job ran the unchanged published input and a second input changing only
Simulink `ode45` `MaxStep` from 0.05 to 0.005 s. The fixture retains time,
heave, heave speed, restoring and excitation force, PTO force, wave elevation,
and the five clock settings for each arm. The source checkout was WEC-Sim
`0753b2e47f2457c078751dcfe5d251d1767b80ab`; Applications was
`d53d4d4c9eda2581f04204f5d394a6ef84bb099e`.

The jointly refined 0.005 s reference is the separate
`../ellipsoid_ode45_maxstep_only/refined_005_common_grid.npz` fixture.
The new audit regenerates its fine-clock arm and holds `dt=MaxStep=0.005`
and `dtOut=0.05` s constant while restoring nonlinear-force and, for CIC,
convolution sampling one clock at a time. These are derived source cases;
the published 0.05 s trajectories remain qualified.

# RM3 PTO Extension fixed-step-only source

`source.npz` contains the 301 samples at the published 0.1 s output grid
from two local MATLAB R2025b Update 7 runs of each pinned
`RM3_PTO_Extension/float` and `/spar` input. The `published` arm leaves
`ode4 FixedStep=simu.dt=0.1` s; the `fixedstep0125` arm changes only
Simulink `FixedStep` to 0.0125 s after `initializeWecSim`. Both retain
`simu.dt=dtOut=nonlinearDt=cicDt=0.1` s, the 30 s no-wave source,
the +5 m PTO extension, and zero PTO spring and damping. Each array holds
time, float heave/speed, spar heave/speed, PTO stroke/speed, internal
mechanics force, and power; the matching `_settings` array stores the five
time settings. The unchanged source outputs reproduce the older published
fixture `../rm3_pto_extension_step.npz` within `1e-7` in logged units.

Pinned WEC-Sim core: `0753b2e47f2457c078751dcfe5d251d1767b80ab`.
Pinned Applications: `d53d4d4c9eda2581f04204f5d394a6ef84bb099e`.
The ordinary Python suite tests this compact source trace, while the
dedicated MATLAB workflow regenerates both arms and applies the same gate
to fresh output. The fixed-step arm is a derived source case; it does not
replace either published 0.1 s trajectory.

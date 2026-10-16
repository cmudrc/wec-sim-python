# RM3 MoorDyn acceleration-delay prefix

`signals.npz` contains six-DOF acceleration signals for both hydrodynamic
bodies in a 2 s MATLAB R2025b Update 7 run. The source was pinned to
WEC-Sim `0753b2e47f2457c078751dcfe5d251d1767b80ab`, Applications
`d53d4d4c9eda2581f04204f5d394a6ef84bb099e`, and the legacy MoorDyn
library `34c0703c9a5b74fd2cc5bb52caaadaaa356d675e` on arm64 macOS.

The temporary case copied published `RM3_MoorDyn_Viz` and changed only
`endTime` from 80 to 2 s, `dtOut` from 0.1 to 0.01 s, ParaView recording
from on to off, `waves.phaseSeed` to 1, and execution mode from accelerator
to normal. It kept the published `simu.dt=0.01` s, ode45 solver, MoorDyn
input, and hydrodynamics. Normal mode avoided the local accelerator SDK
configuration; its body positions matched the archived
[Linux dense source prefix](https://github.com/cmudrc/wec-sim-python/actions/runs/38085085794)
within `1e-12` m. The sea phases and elevation matched exactly.

Simulink port logging captured the acceleration entering and leaving the
hydrodynamic body's `10e-8` s Transport Delay, plus the acceleration entering
its radiation-force block. The 207 logged ode45 major-step times include six
additional times outside the 201 output samples. `output_acceleration` comes
from a separate uninstrumented run of the same case; it equals the logged
delay input exactly at every output time. The signal-logging run used
`ReturnWorkspaceOutputs` to capture `logsout`, so its postprocessing was not
used for this fixture.

The [WEC-Sim developer manual](https://wec-sim.github.io/WEC-Sim/dev/developer/advanced_features.html#added-mass)
describes this tiny delay as an extrapolation of prior acceleration to break
the added-mass algebraic loop.

The portable test checks the active signal connection and Simulink's linear
extrapolation from the two prior major steps. This is a source-timing
diagnostic, not an independent Python motion gate for the published 80 s case.

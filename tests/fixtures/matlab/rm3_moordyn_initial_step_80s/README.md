# RM3 MoorDyn full-duration initial-step audit

`signals.npz` contains both 80 s MATLAB source arms at the published 0.1 s
output interval and all 8,001 native MoorDyn fairlead samples. The only
solver change is ode45 `InitialStep`: `auto` versus `0.01` s. The pinned
application retains `simu.dt=0.01` s, `MaxStep=0.01` s, its 1,000-bin
JONSWAP sea, hydrodynamics, PTO, and MoorDyn. Both arms use phase seed 1,
normal Simulink mode, and disabled ParaView output. These last three settings
match the prior controlled source runs; the unchanged local arm reproduces the
archived published 80 s body positions within `3e-12` m or rad and MoorDyn
loads within `1e-5` N or N m at every common sample.

The source is WEC-Sim `0753b2e47f2457c078751dcfe5d251d1767b80ab`,
Applications `d53d4d4c9eda2581f04204f5d394a6ef84bb099e`, and legacy
MoorDyn `34c0703c9a5b74fd2cc5bb52caaadaaa356d675e`. MATLAB R2025b
Update 7 ran on arm64 macOS with BEMIO-generated RM3 HDF5, SHA-256
`1caecb3061697edbd6c1949092f5b1b1a124d28e1836bc173b226e69e8ba1d00`.
The reproducer is
[`rm3_moordyn_initial_step_80s_audit.m`](../../../matlab/rm3_moordyn_initial_step_80s_audit.m).
The fixture's SHA-256 is
`276a4dcf778e038824099406b44d5c8b75f49511a8c981582056e9fb3bd16cff`.

Changing `InitialStep` alone moves the 80 s MATLAB float/spar surge by
`209.75`/`186.15` mm, shared pitch by `0.0017973` rad, and MoorDyn pitch
moment by `104.05` kN m. Relative to the previously archived fully refined
source, the changed arm has maximum float/spar surge differences of
`6.32`/`7.29` mm, pitch difference `4.77e-5` rad, and mooring pitch-moment
difference `9.64` kN m. The squared-gap reductions are `99.921%`, `99.922%`,
`99.936%`, and `99.611%` respectively. Spar heave has only `84.23%` reduction.

A separate public Python `WEC.run` run advanced motion with the pinned HDF5,
phase values, and live native MoorDyn, without saved source motion or force as
inputs. Compared with the changed 80 s source arm, its maximum float/spar
surge differences are `7.07`/`8.08` mm; shared pitch is `5.60e-5` rad;
MoorDyn surge/heave force and pitch moment are within `1.24`/`0.689` kN and
`8.85` kN m; PTO force is within `706.5` N. The three native fairlead
tensions differ by at most `2.065`, `0.706`, and `0.706` kN. The optional
native test reruns this comparison when `WEC_SIM_RM3_H5`,
`WEC_SIM_MOORDYN_LIBRARY`, and `WEC_SIM_RM3_MOORDYN_LINES` are available.
The source-only test runs in ordinary Python CI.

This attributes a large part of the published source trajectory to ode45
startup-step sensitivity. It is not an independent parity gate for the
unchanged published input or evidence of a WEC-Sim physical-model bug. The
RM3 MoorDyn visualization case remains qualified and its badge count stays
unchanged.

# RM3 MoorDyn initial-step diagnostic

`signals.npz` preserves a controlled two-second comparison of the pinned
`Paraview_Visualization/RM3_MoorDyn_Viz` application. The WEC-Sim source is
`0753b2e47f2457c078751dcfe5d251d1767b80ab`, Applications is
`d53d4d4c9eda2581f04204f5d394a6ef84bb099e`, and native MoorDyn is
`34c0703c9a5b74fd2cc5bb52caaadaaa356d675e`. MATLAB R2025b Update 7
ran both arms on arm64 macOS with the same locally generated `rm3.h5`.

Both arms retain the published `simu.dt=0.01`, 1,000-bin JONSWAP sea,
hydrodynamics, PTO, and MoorDyn. For this bounded diagnostic, they use
`endTime=2`, `dtOut=0.01`, phase seed 1, normal Simulink mode, and no
ParaView output. One arm leaves `ode45 InitialStep=auto`; the other changes
only `InitialStep` to `0.01` after WEC-Sim initialization. Both retain
`MaxStep=0.01`. The reproducer is
[`tests/matlab/rm3_moordyn_initial_step_audit.m`](../../../matlab/rm3_moordyn_initial_step_audit.m).

The saved default arm reproduces the earlier hosted 80-second published-step
trace on all shared times through 2 seconds: identical phase values, wave
elevation within `1e-15` m, and body positions within `1e-12` m or rad.
The two local source arms have exactly identical phases and wave samples.
The checked-in MATLAB reproducer regenerated both source arms with exact
agreement to every saved phase, wave, and body-position value in this fixture.
Changing only `InitialStep` moves float/spar surge by at most 36.419/62.001 mm
and shared pitch by 0.0012433 rad over the two-second window.

An independent public `WEC.run` solve used only the pinned HDF5, source phase
values, case settings, and live native MoorDyn. Against the default MATLAB arm,
its maximum float/spar surge differences were 34.468/60.556 mm and shared
pitch 0.0012679 rad. Against the `InitialStep=0.01` arm, those differences fell
to 1.950/1.446 mm and 0.00002456 rad. The optional native test in
`tests/test_rm3_moordyn_initial_step_portable.py` reruns this independent
comparison when `WEC_SIM_RM3_H5`, `WEC_SIM_MOORDYN_LIBRARY`, and
`WEC_SIM_RM3_MOORDYN_LINES` are supplied. Its source-only test runs in
ordinary Python CI.

This is a short, derived numerical diagnostic. It does not establish
independent parity with the unchanged published 80-second input, and it does
not justify changing the physical Python default or the RM3 badge count.

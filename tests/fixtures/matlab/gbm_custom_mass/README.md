# Barge still-water custom-mass check

`heavy_heave.csv` is the local MATLAB R2025b Update 7 response from the pinned
WEC-Sim core `0753b2e47f2457c078751dcfe5d251d1767b80ab` and Applications
`d53d4d4c9eda2581f04204f5d394a6ef84bb099e`. The starting model is
`Generalized_Body_Modes/GBM.slx` with its published `wecSimInputFile.m` and
`geometry/barge.stl`. It uses the checked-in `../gbm_barge.h5` (SHA-256
`0d98003d54f007f1c50966bb10df723f5bbbf22aa70579360e4ef991a2ceff68`).

For this bounded still-water diagnostic, the only input changes were
`simu.endTime = 1`, `simu.rampTime = 0`, `waves.height = 0`, and
`body(1).mass = 4100000`. The source retained its `ode45` solver and 0.05 s
output step. After adding the pinned core `source` directory to the MATLAB
path, run `wecSim` from that model directory and export:

```matlab
values = [output.bodies(1).time(:), output.bodies(1).position(:,3), ...
          output.bodies(1).velocity(:,3)];
writematrix(values, 'heavy_heave.csv');
```

The 21 CSV rows are time, heave, and heave speed. File SHA-256 is
`7571786ae213a63870f0f8f1580b642ccaf1b586c42c59c172bb298c7db3caca`.
The independent public Python run differs by at most 2.85 mm in heave and
7.99 mm/s in heave speed across this 1 s window. This diagnoses the static
weight/buoyancy term; it does not alter the published-case badge count.

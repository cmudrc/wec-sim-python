# Derived RM3 solver-step negative control

This fixture is a **10 s derived experiment**, not a published WEC-Sim input.
It uses the pinned core `examples/RM3/RM3.slx`, which has no MoorDyn block.
The copied input sets the visualization application's 2 m, 8 s seeded
Traditional JONSWAP sea, 70 m water depth, zero ramp, initial spar heave of
−0.21 m, 1.2 MN s/m PTO damping, and `ode45`. Both source runs keep
`simu.dt=0.01` s; only Simulink `MaxStep` changes from `0.01` to
`0.0003125` s. They record 1,001 common samples. The plotting callback is
disabled in the disposable copy after simulation. The pinned checkout stays
untouched.

The two MATLAB runs have **identical** 1,000 phase values and all 1,001 wave
samples. Refining `MaxStep` changes float and spar surge by up to 265.844
and 239.644 mm, and their common pitch by up to 0.0016921 rad. Thus a large
coarse-step RM3 surge shift is possible without MoorDyn. This does not prove
that every part of the published MoorDyn visualization gap has the same cause:
the core RM3 Simulink model differs from `RM3MoorDyn.slx`.

The portable test independently advances the public Python two-body RM3 model
using only the pinned HDF5 and saved sea phases. Against the fine-step MATLAB
run, maximum float/spar surge differences are 0.244/0.132 mm, and all active
translational speeds differ by less than 0.090 mm/s. It reuses the existing
RM3 MaxStep-only gate limits: 2 mm translational position, `5e-5` rad pitch,
0.7 mm/s translational speed, and `5e-5` rad/s pitch speed. The published
MoorDyn case remains qualified; no case badge changes.

## Provenance and reproduction

- WEC-Sim core: `0753b2e47f2457c078751dcfe5d251d1767b80ab`.
- MATLAB: R2025b on mali, normal Simulink mode.
- Original `examples/RM3/wecSimInputFile.m` SHA-256:
  `ce044db7695c129d352920c61ec45371860e03c70df2b19bf72097e4e26a6959`.
- `examples/RM3/RM3.slx` SHA-256:
  `9aab8875d9f35b408b3cca1eabdd95e4219b2e6f94dddf19c172fba729576c63`.
- Both generated `wecSimInputFile.m` files SHA-256:
  `1cff13fe8f62c9e391bbdeb56e11cb690d66ccb992afe37a0ed669bc81a7972e`.
- [`rm3.h5`](rm3.h5) is byte-identical to the pinned core input, SHA-256:
  `a909360c6479773c113f792be72afb1351cae9299d62e501fa80561dca450aa3`.
- [`source.npz`](source.npz) stores the coarse and fine phase, wave, and both
  body histories, SHA-256:
  `2ce203b8c6a6b90d1d2439618f71ab71b57eda3eba13e6b04afdff4998ea6d02`.

From this repository checkout, invoke
[`rm3_unmoored_step_audit.m`](../../../matlab/rm3_unmoored_step_audit.m) in MATLAB's
base workspace with `rm3CoreRoot` set to the pinned WEC-Sim checkout and
`rm3OutputRoot` set to a **new, nonexistent** output directory. For example:

```matlab
rm3CoreRoot = '/path/to/WEC-Sim';
rm3OutputRoot = '/tmp/rm3-unmoored-audit-new';
run('tests/matlab/rm3_unmoored_step_audit.m');
```

The script writes `coarse_*.csv` and `fine_*.csv`. `source.npz` is a compressed
bundle of those eight CSV arrays, loaded with `numpy.loadtxt` and saved under
their filename stems with `numpy.savez_compressed`. The paired Python check is
`tests/test_rm3_unmoored_step_diagnostic.py`.

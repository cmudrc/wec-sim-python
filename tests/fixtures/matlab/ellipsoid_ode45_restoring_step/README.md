# Ellipsoid copied-model continuous nonlinear MaxStep continuation

These fixtures preserve the 150 s Regular and RegularCIC source traces from
[GitHub Actions run 40131187182](https://github.com/cmudrc/wec-sim-python/actions/runs/40131187182),
which succeeded at Python-repository commit
`b910e47d5a6ad2edd5b6af183374d83d74db1e8a`. The workflow pinned WEC-Sim
core `0753b2e47f2457c078751dcfe5d251d1767b80ab` and Applications
`d53d4d4c9eda2581f04204f5d394a6ef84bb099e`, using MATLAB R2025b. The
[runner](../../../matlab/ellipsoid_ode45_restoring_timing_smoke.m) copied each
published application input and model, made the restoring, Froude–Krylov, and
wave-elevation blocks continuous, then varied only the copied model's `ode45`
MaxStep among these arms. The published input stays at 0.05 s. The physical, output,
nonlinear, and CIC clocks stay at 0.05 s, with the original 150 s duration,
50 s ramp, 4 m/6 s sea, and 1.2 MN s/m PTO damper.

Each NPZ contains 3,001 output samples for four copied-model MaxSteps:
0.025, 0.0125, 0.00625, and 0.003125 s. It stores heave, speed, restoring
force, PTO stroke/speed/force/absorbed power, incident wave, nine source
settings, and the common time grid. The separate
[full force-timing fixture](../ellipsoid_ode45_restoring_timing/README.md)
supplies the 0.05 s continuous nonlinear base. The packer checked all eight
fresh source HDF5 files against the 38 datasets in the archived source HDF5
by value, including NaN locations. Binary HDF5 hashes need not match.

The original predeclared screen applied to the 0.025→0.0125 s change. Its
pointwise limits were 1.0 mm heave, 1.3 mm/s speed, 1.6 kN PTO force, and
1.3 kW power for Regular; 1.1 mm, 1.5 mm/s, 2.0 kN, and 1.6 kW for
RegularCIC. Both had a 0.1% signed-work change limit and required contraction
by a factor below 0.7 from the preceding change. Every measure contracted,
but the **predeclared screen failed**: Regular exceeded speed, force, power,
and work limits; RegularCIC exceeded all five limits.

The following 0.00625 and 0.003125 s runs are a **post-screen diagnostic**.
Their last changes fall below the old numerical ceilings, but this does not
retroactively pass the predeclared 0.0125 s screen.

| Case | MaxStep pair (s) | Heave (mm) | Speed (mm/s) | PTO force (kN) | Power (kW) | Signed-work change |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Regular | 0.05→0.025 | 1.700 | 3.232 | 3.879 | 3.287 | 0.4186% |
| Regular | 0.025→0.0125 | 0.907 | 1.722 | 2.066 | 1.339 | 0.2238% |
| Regular | 0.0125→0.00625 | 0.440 | 0.797 | 0.956 | 0.693 | 0.1125% |
| Regular | 0.00625→0.003125 | 0.242 | 0.395 | 0.474 | 0.367 | 0.0634% |
| RegularCIC | 0.05→0.025 | 2.218 | 4.511 | 5.413 | 4.598 | 0.5494% |
| RegularCIC | 0.025→0.0125 | 1.115 | 2.356 | 2.828 | 1.940 | 0.2785% |
| RegularCIC | 0.0125→0.00625 | 0.556 | 1.082 | 1.299 | 0.968 | 0.1410% |
| RegularCIC | 0.00625→0.003125 | 0.293 | 0.549 | 0.659 | 0.504 | 0.0712% |

An independent Python run uses the archived source HDF5/STL and configured
wave/PTO input, without source motion or force histories. At the diagnostic
0.003125 s arm, all 3,001 pointwise samples and signed absorbed work fall
within the inherited Python physical limits:

| Case | Heave gap (mm) | Speed gap (mm/s) | PTO force gap (kN) | Power gap (kW) | Signed-work gap |
| --- | ---: | ---: | ---: | ---: | ---: |
| Regular | 0.220 | 0.383 | 0.460 | 0.371 | 0.0593% |
| RegularCIC | 0.406 | 0.516 | 0.619 | 0.412 | 0.0700% |

The inherited physical limits are 3 mm/4 mm/s/5 kN/4 kW for Regular and
3.5 mm/4.5 mm/s/6 kN/5 kW for RegularCIC, plus 0.1% signed work for both.
This result concerns only the copied continuous nonlinear model. The
unchanged published-input gate remains red and its status is unchanged.

`regular.npz` is 533,989 bytes (SHA-256
`b19a80a8d2451689056d2dea4ba85ecbaa0961e9cc850a641c4c6271df589246`);
`regular_cic.npz` is 633,842 bytes (SHA-256
`647c1be3f837912442e380417a98af00245bb7ba1b0da3f9d6c7348d1b3e1546`).
Regenerate and compare the archives from the workflow artifacts:

```sh
gh run download 40131187182 -D .evidence/force_timing_finer
.venv/bin/python tools/pack_ellipsoid_ode45_restoring_step.py Regular .evidence/force_timing_finer .evidence/regular-step-repacked.npz --expect tests/fixtures/matlab/ellipsoid_ode45_restoring_step/regular.npz
.venv/bin/python tools/pack_ellipsoid_ode45_restoring_step.py RegularCIC .evidence/force_timing_finer .evidence/regular-cic-step-repacked.npz --expect tests/fixtures/matlab/ellipsoid_ode45_restoring_step/regular_cic.npz
.venv/bin/python -m pytest -q -s tests/test_ellipsoid_ode45_restoring_step_portable.py
```

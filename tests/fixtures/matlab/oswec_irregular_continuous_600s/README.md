# OSWEC irregular continuous-heading control

This is a **derived 600 s control**, not the published `Variable_Hydro/Passive_Yaw`
selected-bank input. The published irregular case remains qualified, and this
control does not increase the OSWEC badge count.

The direct MATLAB runs use WEC-Sim 7.1.0 source commit
`0753b2e47f2457c078751dcfe5d251d1767b80ab` and Applications commit
`d53d4d4c9eda2581f04204f5d394a6ef84bb099e` on local R2025b. The
Applications `Variable_Hydro/Passive_Yaw/OSWEC.slx` model is byte unchanged
(SHA-256 `80eb2552f2cd36c96c6dc362d8a08878e9ada505a8b8437feaa574d616642099`).
The source interpolation function `irregExcFullDirF.m` is also unchanged in the
baseline runs (SHA-256 `9b16fe36ba655b700662765655ea6411f795f1b8bff4a847d318c976448cadba`).

The archived input files change only the flap hydrodynamic source from a bank
of single-direction HDF5 files to the original all-direction OSWEC HDF5,
select WEC-Sim's existing continuous passive-yaw interpolation with threshold
zero, and use 0.01 or 0.005 s fixed `ode4` steps. Geometry, mass, inertia,
wave, seed, ramp, 600 s duration, radiation memory, and PTO remain as in the
published application. `calcIndex.m` supplies a constant unused index because
the copied Simulink model still contains the indexing block. The application
plotting hook is disabled because it assumes a selected-bank index. No source
states, forces, or headings are fed into the Python run.

The local MATLAB HDF5 SHA-256 is
`a423b49b9e09e3303902e576a1f57a696f01c51e877d3271a7d71b59e00493e5`.
The corresponding [Python fixture](../oswec_multi_wave/oswec.h5) stores the
same numeric hydrodynamic arrays; the largest dataset difference is
`1.46e-11` in an IRF array from HDF5 regeneration. `gate_spec.json` records
the unchanged issue #621 yaw, speed, torque, and PTO-work limits before either
full-duration direct result was inspected.

The 40 s calibration run at 0.005 s matched independent physical `WEC.run`
within `5.84e-8` rad yaw, `3.66e-8` rad/s speed, `0.00439` N m torque, and
`0.000627%` PTO work. A separate MATLAB accelerator run produced identical
saved values but was only 3% faster than normal mode, so the 600 s baselines
use normal mode. `source.npz` stores the two direct full-duration runs after
they complete; the portable test gates MATLAB step change and the public
physical Python trajectory separately.

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

The pinned Applications checkout contains WAMIT data and BEMIO scripts, not
the generated HDF5 files. Generate the full-direction `oswec.h5` with
`_Common_Input_Files/OSWEC/hydroData/bemio.m` and the single-direction
`oswec_10.h5` with `Variable_Hydro/Passive_Yaw/hydroData/bemio.m`. In that
script's final writer loop, write only the entry whose `file` ends in
`_10`; the other selected headings are unused by this control. To reproduce
either source trace, copy the pinned application's `Variable_Hydro/Passive_Yaw`
directory next to its `_Common_Input_Files` directory. Replace
`wecSimInputFile.m` with the
corresponding archived step file and copy `calcIndex.m` and
`userDefinedFunctions.m` into that case.
Call `run_source(source_dir, case_dir, output_dir, step)` from MATLAB. It
exports time, flap yaw and speed, signed PTO torque and power, and incident
wave elevation without modifying the WEC-Sim or Applications checkouts.

The local MATLAB HDF5 SHA-256 is
`a423b49b9e09e3303902e576a1f57a696f01c51e877d3271a7d71b59e00493e5`.
All 34 common non-excitation numeric datasets in that full-direction file and
the published `oswec_10.h5` bank member are exactly equal, including added
mass, radiation, restoring stiffness, and body properties.
The flap excitation real and imaginary arrays also match exactly at the
initial 10° incidence; only the heading-selection law changes during yaw.
The corresponding [Python fixture](../oswec_multi_wave/oswec.h5) has all 68
dataset paths and the same NaN masks. Only four regenerated IRF arrays differ;
their largest absolute difference is `1.46e-11`. `gate_spec.json` records
the unchanged issue #621 yaw, speed, torque, and PTO-work limits before either
full-duration direct result was inspected.

The 40 s calibration run at 0.005 s matched independent physical `WEC.run`
within `5.84e-8` rad yaw, `3.66e-8` rad/s speed, `0.00439` N m torque, and
`0.000627%` PTO work. A separate MATLAB accelerator run produced identical
saved values but was only 3% faster than normal mode, so the 600 s baselines
use normal mode. The two direct 600 s source runs took 3,658 and 6,997 s.
`source.npz` stores both full-duration traces (SHA-256
`a282fd2726727917ccd3de19f3ed683536541c5504055d7428d2c6fe7554dea8`);
`evidence.json` records the measured gaps and input hashes. The portable test
gates MATLAB's own step change and the public physical Python trajectory
separately.

MATLAB 0.01-to-0.005 s self-change is `3.21e-5` rad yaw,
`4.85e-6` rad/s speed, `0.582` N m torque, and `0.0121%` PTO work. Against
direct 0.005 s MATLAB, independently advanced physical Python differs by
`1.25e-5` rad yaw, `1.92e-6` rad/s speed, `0.231` N m torque, and
`0.00492%` work; its incident wave differs by `2.13e-13` m. Both comparisons
pass the unchanged issue #621 limits over all 600 s. All three independently
computed relative headings remain between `-8.756°` and `10.270°`, inside
the original `-30°` to `30°` selected-bank span. The MATLAB work is
2,211.017 J at 0.01 s and 2,211.284 J at 0.005 s; Python gives 2,211.175 J.
This derived control does not promote the published nearest-bank case.

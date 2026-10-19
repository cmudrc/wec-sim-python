# Sphere pressure-integration mean drift

This is the 30 s PM-wave Sphere variant in
[sphere_irregular_mean_drift](../sphere_irregular_mean_drift/README.md), with
only `body(1).meanDrift` changed from control surface (`1`) to pressure
integration (`3`). The pinned HDF5 contains both coefficient families. A local
R2025b run used core `0753b2e47f2457c078751dcfe5d251d1767b80ab` and
Applications `d53d4d4c9eda2581f04204f5d394a6ef84bb099e`. The shared
[MATLAB generator](../../../matlab/sphere_irregular_mean_drift_reference.m)
reproduced the independently prepared source run exactly across all body,
wave, and 500 component samples.

`source.npz` stores all 3,001 body and wave samples, sea components, and body
mass properties. Its SHA-256 is
`ac602a48bd02e017b2fdc4501ae7f8911d7f8afd970c64e1b3460f7f547fd608`.
The independent public Python `WEC.run` has maximum MATLAB differences of
`5.1e-15` m wave elevation, `1.3e-10` N or N m excitation, 0.0425/0.134 mm
surge/heave, `5.2e-8` rad pitch, and 0.425 mm/s speed. Against the same-sea
control-surface MATLAB run, pressure integration shifts final surge by
18.95 mm and pitch by 0.00329 rad. This gates a distinct selectable coefficient
path; it is a derived input and does not add a published-case badge.

# Sphere momentum-conservation mean drift

This is a derived variant of the pinned Applications `Mean_Drift` input. The
only input edit is `body(1).meanDrift = 2` in place of `1`; wave, geometry,
body properties, solver, time step, and 100 s duration are unchanged. A local
MATLAB R2025b Update 7 run on 2026-10-15 used WEC-Sim core
`0753b2e47f2457c078751dcfe5d251d1767b80ab` and Applications
`d53d4d4c9eda2581f04204f5d394a6ef84bb099e`. Pinned BEMIO regenerated
`sphere.h5` from the Applications WAMIT files.

`source.npz` retains the full 10,001-sample MATLAB body and wave exports,
the six-component real, imaginary, and mean-drift excitation coefficients,
and the four mass properties. Its SHA-256 is
`bc1a32a4b1e5f78b35396952754b80517825959553aff9898279adf9a76e768c`.
The regenerated HDF5 frequency grid, directions, displacement, restoring,
excitation, infinite-frequency added mass, and momentum-drift coefficients
match the existing [Sphere mean-drift HDF5](../sphere_mean_drift/sphere.h5)
exactly; the active radiation IRF differs by at most `1.07e-14`. The fitted
state-space realization differs but is not used by this convolution case.

The public `WEC.run` gate independently advances the full 100 s from that
existing HDF5. It reuses the published control-surface case's motion and
force limits. Observed maximum MATLAB/Python differences were `4.21e-5` m
surge, `5.51e-4` m heave, `4.62e-8` rad pitch, `8.05e-5` m/s surge speed,
`1.74e-3` m/s heave speed, `1.32e-7` rad/s pitch speed, `2.85e-15` m wave
elevation, and `9.30e-12` N or N m excitation. The mean-drift coefficient
force identity agrees within `3.56e-15` N.

This checks a distinct physical coefficient path. It is not another
published input, so the Sphere case badge remains unchanged.

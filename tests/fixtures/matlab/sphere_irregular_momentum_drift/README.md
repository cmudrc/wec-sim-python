# Sphere momentum-conservation drift in a PM sea

This is a 30 s derived variant of the pinned Applications `Mean_Drift` input.
The [MATLAB generator](../../../matlab/sphere_irregular_mean_drift_reference.m)
uses the same 0.5 m, 8 s PM sea, seed 1, HDF5, three-coordinate floating
joint, and 0.01 s step as the control-surface and pressure-integration
variants, selecting only `body(1).meanDrift = 2`. Local MATLAB R2025b used
WEC-Sim core `0753b2e47f2457c078751dcfe5d251d1767b80ab` and
Applications `d53d4d4c9eda2581f04204f5d394a6ef84bb099e`.

`source.npz` stores the 3,001-sample body and wave outputs, all 500 spectral
components, and body mass/inertia as float64 without rounding. Its SHA-256 is
`9faa78a7bcc723a78eb0c210fbf71666600254fda352c6f60f2edccfd100ca6c`.
The source CSV hashes are `50ce8f3849384be4429698b99574879464056fffbb776d7708917c336d684699`
for body, `5bcff824f98052e14e4d7e804fc826fa49eaff470fe033c2f8f1daf6b7658bf1`
for wave, `f0b27cc134cd3b55a19485289e9335b7c9ecfa666f97ef0d17d296b9f40dae53`
for components, and `6d9b981cd5779dd436c9df42d0ca81d9109f5ff3b08b4003791fceedc4225921`
for mass. The generated input SHA is
`b7f6cbf5d9c37edc32b07c5545e069c1fe58f68f2ecc174ba88f9205cb21b3aa`.
An independent second local R2025b run reproduced all four CSV hashes exactly.
The source-matched HDF5 is the existing
[`sphere_mean_drift/sphere.h5`](../sphere_mean_drift/sphere.h5), SHA-256
`418f7fa2059ebe1d00f8520c55962dbafc6f256a8d26f4117afcc3d495abe3fb`.

The independent public `WEC.run` test regenerates every PM frequency, width,
amplitude, and phase. Maximum MATLAB/Python gaps are `5.05e-15` m wave
elevation, `1.30e-10` N or N m across six excitation channels, 0.0425 mm
surge, 0.134 mm heave, `5.16e-8` rad pitch, and 0.425 mm/s active speed.
An HDF5 spectral-energy calculation reproduces the reported drift force
within `1.2e-16` N. The original control-surface and new momentum source
traces end at 0.025389 and 0.000197 rad pitch, respectively, confirming that
the coefficient choice changes the simulated response. Existing
control-surface motion and force limits remain unchanged; this test adds a
separate nonzero momentum-surge check and verifies its zero non-surge drift
coefficients.

This derived configuration is not another published input. The Sphere badge
and 63-case denominator are unchanged.

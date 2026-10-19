# Sphere spectral mean drift

This is a derived 30 s PM-wave variant of the published Applications
`Mean_Drift` case. It keeps the Sphere model, control-surface mean-drift
selection, body properties, 0.01 s step, and 10 s radiation memory. The input
changes wave height to 0.5 m, peak period to 8 s, phase seed to 1, ramp to
10 s, and end time to 30 s. The source's optional regular-wave plotting script
is removed from the copied case because it assumes a scalar amplitude.

The local R2025b source run used pinned WEC-Sim core
`0753b2e47f2457c078751dcfe5d251d1767b80ab`, Applications
`d53d4d4c9eda2581f04204f5d394a6ef84bb099e`, and the previously pinned
[`sphere.h5`](../sphere_mean_drift/sphere.h5) from that application's BEMIO
input. [The source script](../../../matlab/sphere_irregular_mean_drift_reference.m)
exports all 3,001 body, wave, and excitation samples plus the 500 spectral
components. `source.npz` is SHA-256
`e620338b9f34f5ad44f3b2012d6a6a7e1397c8a9992e2c692218d30c6a2ad475`.

The public Python `WEC.run` independently reproduces the source sea within
`5.1e-15` m, all six excitation components within `1.4e-10` N or N m, and
surge, heave, and pitch within `0.043` mm, `0.134` mm, and `5.2e-8` rad over
30 s. The largest speed difference is `0.425` mm/s. Removing mean drift from
the same Python run changes final surge by `0.081` m and pitch by `0.0252`
rad, so the case exercises the selected spectral coefficients. This is a
derived feature gate, not a new published input or a change to the 63-case
badge.

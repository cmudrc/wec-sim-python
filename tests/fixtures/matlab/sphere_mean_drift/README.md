# Published Sphere mean drift

`source.npz` preserves the published 100 s `Mean_Drift` MATLAB output from the
successful [R2025b source run](https://github.com/cmudrc/wec-sim-python/actions/runs/38677803196)
(`SPHERE_MEAN_DRIFT`, run commit `48e11f79682346404acb3762fc80d10ddcf8a34c`).
The workflow pins WEC-Sim core
`0753b2e47f2457c078751dcfe5d251d1767b80ab` and Applications
`d53d4d4c9eda2581f04204f5d394a6ef84bb099e`.

The four arrays are loaded as float64 from the source CSVs and stored without
rounding or interpolation using `numpy.savez_compressed`:

| Archive key | Source CSV suffix | SHA-256 |
| --- | --- | --- |
| `body` | `Mean_Drift_body1` | `f6f7fef648a37798bd7e3b65589eed501ca96234533f6d3e6d64f92f6e5c4a79` |
| `wave` | `wave` | `ed5d59834453340b2e3b6e2810431bddab961d9212578f3deed69ceb4d85e489` |
| `coefficients` | `excitation_coefficients` | `7a4612385468616cb35313f36bfa9b271b0551abee002dff8c677ee7adc8c4a8` |
| `mass` | `mass_properties` | `6d9b981cd5779dd436c9df42d0ca81d9109f5ff3b08b4003791fceedc4225921` |

The archive SHA-256 is
`b6d8bfaf39a898d18db7274f09120bafdf3bd9d4ede941b4de0195a8d52128c0`.
`sphere.h5` is a byte copy of that run's
`applications/Mean_Drift/hydroData/sphere.h5` (SHA-256
`418f7fa2059ebe1d00f8520c55962dbafc6f256a8d26f4117afcc3d495abe3fb`).
It differs from `examples/data/sphere.h5`, so the source-matched file is
required for this comparison.

Ordinary Python CI runs the existing full-duration public `WEC` motion,
excitation, drift, and radiation-force gates with this fixture. The live MATLAB
workflow still uses freshly generated CSVs and its application HDF5. This
makes an existing published-case gate portable; it does not add a case badge.

# Sphere variable-mass source fixture

`source.zip` contains the nine original `draft1.h5`–`draft9.h5` files and five
unmodified CSV exports from the successful [R2025b 900 s source run](https://github.com/cmudrc/wec-sim-python/actions/runs/38282150325)
(`matlab-SPHERE_VARIABLE_MASS-baseline`, repository commit
`7c403d3e2e07b9940e4a2a280d375000b8985fea`). That workflow pins WEC-Sim
core `0753b2e47f2457c078751dcfe5d251d1767b80ab` and Applications
`d53d4d4c9eda2581f04204f5d394a6ef84bb099e`. The ZIP SHA-256 is
`9dc96155bad3f2b70155839d263cefa463cd28dd4134f220e830daa6139c9715`.

The CSVs are `SPHERE_VARIABLE_MASS_Variable_Mass_body1.csv`,
`SPHERE_VARIABLE_MASS_Variable_Mass_pto1.csv`, `SPHERE_VARIABLE_MASS_active_state.csv`,
`SPHERE_VARIABLE_MASS_states.csv`, and `SPHERE_VARIABLE_MASS_wave.csv`. The ZIP
retains their original artifact paths. Their SHA-256 hashes, in that order,
are:

```text
2cc523581d1e61f37874370e2235858d66bd37c6a2fc3f9f5b1e26b6bed96776
193f103172362f0388f1287fa27182a02ad4f10d9786674ff78c86d8638ad0b6
1d2bf089bd3f9921a05f7f1e115776257aac7331cdda5a96ad1ba918f1d5213d
fec34311189731cb148ae1f748623344bdf953a51d582bee8d2e148b2015ec93
ef476d73829f53bba92940d34f9f56cc88cd0bf05126711e008a5ea77a63d31e
```

The `draft1.h5` through `draft9.h5` SHA-256 hashes, in order, are:

```text
10e020bc87c49a4cc8303ac1eb67bacf38b403ddd4c4dcd5a76f7e4c578ec7d0
24614a323670d9a2f28c20562d51d1f87ec9293e86cf89191d955efaa673a778
1e974359fd2492db8df1231acd3632ad94b479d8c728ca8a482da03a4869b531
839b4cd2559a0a943437fb501c50d90a3bee437c26a3432e190b69c45ff1d885
52d914500142c419b23fb12e3f2ad8d46571dbe5a73c546dfa27b8779c99d70f
3a02e04a28e57d87c41996d520c644baaea62ee3209b47edd76fb693ead1a7c9
e20b75b6d1da9437ec37ba3d1aecc8bf043893709ccd099014b9e141d3c8d73c
3d639ff556988da4b2db91d86c4d1fa6e0b846466aedcb8e72ba25bd79f3fdbe
0f7dbd7fc43a999e24db95132d61445bf5fe42126ff1368815d2101a50b736d5
```

The ordinary Python test extracts these files into a temporary directory,
checks source state selection and force assembly at all 9,001 source output
times, then independently advances the public `WEC.variable_body` for the full
900 s. It uses the same motion and force limits as the fresh MATLAB workflow.
The additional PTO speed, pointwise absorbed-power, and integrated-work gates
are 0.005 m/s, 2.5 kW, and 1%, respectively; the measured maxima are
0.004135 m/s, 1.052 kW, and 0.00376% work. This tightens an already-gated
case; it does not add a new badge case.

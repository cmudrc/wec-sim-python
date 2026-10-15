# Sphere sampled-elevation reference

This is a derived 40 s input from the pinned Sphere `Free_Decay/0m` model:
the no-wave input is replaced by a two-frequency `elevationImport` record.
The [R2025b source run](https://github.com/cmudrc/wec-sim-python/actions/runs/37785332545)
used run commit `3b8b5978cb98b18e6b5c6676a35816545da184d1` and the pinned
WEC-Sim core and Applications sources. It is a sampled-sea dynamics check,
not an additional published reference-model case.

`etaData.mat` is the original 4,001-row imported record; its SHA-256 is
`c15c892166f8ac1e50812380838c142700e5768a0fc3b16926a39541879d7214`.
`source.npz` retains every source sample without rounding: body CSV columns
0, 3, 9, and 19:25 provide time, heave, heave speed, and six excitation
components; the wave CSV provides sampled elevation. The source body and
wave CSV SHA-256 values are
`9428d69be24bf0d41cd50ce609e83e73c406d50ca62d9eb4c600145e585d591e`
and `068e0061f06cbdb569f56160e22ea5d5e7f505cde699a2ddf504455cf92a5c9d`.
The compact NPZ SHA-256 is
`583e881e6c604ff8c395b85eec7b61deb90ccfd96c106cb74efc66ff1499910c`.
All 46 datasets in the generated Sphere HDF5 match the bundled
`examples/data/sphere.h5`, which the portable test reads. The live-source
path still reads the original HDF5, MAT, and CSV files when both source
environment variables are supplied. Both paths apply the same full-duration
motion, excitation, and elevation gates.

# Published OSWEC two-sea fixture

The successful [R2025b OSWEC_MULTI_WAVE job](https://github.com/cmudrc/wec-sim-python/actions/runs/38130034731)
uploaded `matlab-OSWEC_MULTI_WAVE-baseline`. It used WEC-Sim core
`0753b2e47f2457c078751dcfe5d251d1767b80ab` and Applications
`d53d4d4c9eda2581f04204f5d394a6ef84bb099e`. The job added phase seeds
1 and 2 to the published two-sea input and ran its full 100 s duration.

`source.npz` retains the eight numerical CSVs from that artifact as float64
without rounding: both 500-bin component sets, both 1,001-sample elevation
histories, both body histories, the PTO history, and the source settings.
Its SHA-256 is `43c9afd588beeb1a837f0f670f33e4bf8d88340b75f2db4c905991a8b7abd7e4`.

`oswec.h5` is the artifact's BEMIO HDF5 repacked with gzip and shuffle.
All 68 datasets retain their original shape, dtype, values (including NaNs),
and attributes. The original file's SHA-256 is
`b77d918594655d01fa6bda65a651646b20326c2e9a8a01b30cbda7f314b41604`;
the repacked file's SHA-256 is
`b9ccf3a32e0893c23b1645c15a2baef4a0c6f678b90d177d5797e3c65de4b7d3`.
The bundled `examples/data/oswec.h5` has different coefficients and cannot
replace this published-case input.

`tests/test_oswec_multi_wave_portable.py` independently builds both seas from
the saved phases, advances the public hinge over 100 s, and checks wave,
both-body excitation, motion, PTO torque, and mechanical power against the
existing published-case limits. This makes an already counted OSWEC case
durable in ordinary Python CI; it does not add a case badge.

# MOST blade-load reference

`source.npz` retains the blade-element input arrays, five hub and wind states,
and three-blade root loads from the pinned MATLAB MOST BEM gate. It is a
165,840-byte compressed export of `matlab-most-bem.mat` from [GitHub Actions
run 38264352282](https://github.com/cmudrc/wec-sim-python/actions/runs/38264352282).
That run uses WEC-Sim core `0753b2e47f2457c078751dcfe5d251d1767b80ab`,
Applications `d53d4d4c9eda2581f04204f5d394a6ef84bb099e`, and MATLAB
R2025b. The source data came from `tests/matlab/most_bem_baseline.m`.

The archived MAT file has SHA-256
`5891f30f8cd48cc32417752b4624b581f14508f24537f57f7aa1241e7d853e6f`;
this NPZ has SHA-256
`d44d45392cddb25806033fbd64fafb33e44df5445fd9617374827d3951afd61e`.
The export maps each `bem.bladedata_*` blade geometry field to the matching
`MostBEM` constructor field, subtracts one from MATLAB's one-based airfoil
indices, and retains the source solver parameters. `hub_state`, `wind_speeds`,
`blade_pitch`, and `root_loads` come directly from `q`, `wind_speeds`,
`bladepitch`, and `loads` in the MAT file.

The ordinary-CI gate independently computes all six root loads for all three
blades at each of the five states. The original MATLAB workflow still checks
blade-table parsing and regenerates the source comparison. This fixture does
not establish the published MOST coupled-platform trajectory by itself.

# RM3 seeded two-column spectrum motion

Pinned WEC-Sim `0753b2e47f2457c078751dcfe5d251d1767b80ab` ran its
`examples/RM3/RM3.slx` with the seven-row, two-column
[`spectrum.mat`](../imported_spectrum_two_column/spectrum.mat), phase seed 1,
`ode4` at 0.1 s, a 1 s ramp, and a 30 s end time. The five BEM-range bins
remain after filtering. MATLAB R2025b Update 7 on mali generated the source
motion in `source.npz`; the reproducible script is
[`rm3_two_column_motion_reference.m`](../../../matlab/rm3_two_column_motion_reference.m).
The source model used its bundled `examples/RM3/hydroData/rm3.h5` with SHA-256
`a909360c6479773c113f792be72afb1351cae9299d62e501fa80561dca450aa3`.
Python uses the repository's bundled RM3 HDF5, whose hydrodynamic data give
nearly identical 30 s motion despite a different file hash.

`source.npz` stores the full 301-row wave, two body pose/velocity, PTO, and
five-bin frequency/phase records without rounding. Its SHA-256 is
`4b590fa83fde281abe2ce839409c730aa6bd58b81358fe52cf5477b6e34855e9`.
The public Python gate uses only the spectrum and bundled HDF5 as simulation
inputs. This is a derived 30 s RM3 sea, not another published reference case.

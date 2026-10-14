# RM3 body-to-body published cases

`source.npz` contains the full 4,001-sample body and PTO CSV exports for
all six published RM3 B2B cases. Cases 1 and 2 use regular waves; Cases 3
and 4 use regularCIC waves and convolution radiation; Cases 5 and 6 use
regularCIC waves and fitted state-space radiation. Each pair varies
body-to-body hydrodynamics. They came from the successful
[pinned R2025b baseline run](https://github.com/cmudrc/wec-sim-python/actions/runs/38611588704)
on 2026-10-14, with WEC-Sim core commit
`0753b2e47f2457c078751dcfe5d251d1767b80ab` and Applications commit
`d53d4d4c9eda2581f04204f5d394a6ef84bb099e`. The archive was made by
loading all 18 corresponding source CSVs as float64 and calling
`numpy.savez_compressed` without rounding. Its SHA-256 is
`ca391f8929df2810a2670fdfe4804f0fa58dee8979f71849c26c17ab5818b1e4`.

The portable gates run each published 400 s input from the bundled RM3 HDF5.
Cases 1 and 2 compare regular-wave body and PTO motion. Cases 3 and 4 compare
the physical implicit added-mass scheme and the opt-in
Simulink delay scheme. Cases 5 and 6 compare the opt-in fitted-radiation
solver, and a separate diagnostic verifies that the pinned fit changes the
source motion relative to convolution. These are numerical comparisons: the
pinned state-space fit has negative damping in active low-frequency joint
modes. Body surge, heave, pitch, PTO stroke, speed, force, and power are
checked at every saved time against the existing predeclared limits.
The external Applications/MATLAB-output path remains available for fresh
source runs. This fixture makes all six already gated published cases durable in
ordinary Python CI; it does not add to the case badge count.

The pinned Applications checkout generates `rm3.h5` from `rm3.out` with
`bemio.m`; it does not commit the HDF5. A local R2025b regeneration with the
pinned core reproduced both bodies' restoring, added-mass, excitation, and
frequency-domain damping arrays exactly. The largest radiation-kernel
difference from the bundled HDF5 was `2.73e-12`; displaced volumes differed
by `0.001` and `0.004` m³. Its state-space matrices are different realizations,
but the existing fitted-transfer gate passes. All ten full-duration comparisons
also passed using this regenerated HDF5 with the same archived MATLAB outputs.
Ordinary CI uses the bundled HDF5 to avoid requiring MATLAB at test time.

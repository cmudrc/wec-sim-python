# RM3 body-to-body regularCIC published cases

`source.npz` contains the full 4,001-sample body and PTO CSV exports for
published `B2B_Case3` and `B2B_Case4` (convolution radiation), and
`B2B_Case5` and `B2B_Case6` (fitted state-space radiation). Each pair varies
body-to-body hydrodynamics. They came from the successful
[pinned R2025b baseline run](https://github.com/cmudrc/wec-sim-python/actions/runs/38611588704)
on 2026-10-14, with WEC-Sim core commit
`0753b2e47f2457c078751dcfe5d251d1767b80ab` and Applications commit
`d53d4d4c9eda2581f04204f5d394a6ef84bb099e`. The archive was made by
loading the six corresponding source CSVs as float64 and calling
`numpy.savez_compressed` without rounding. Its SHA-256 is
`6257210c8f7da1e8311ee16d52573d870cc2fceb361850a4966037cb203c7bb0`.

The portable gate runs each published 400 s input from the bundled RM3 HDF5.
Cases 3 and 4 compare the physical implicit added-mass scheme and the opt-in
Simulink delay scheme. Cases 5 and 6 compare the opt-in fitted-radiation
solver, and a separate diagnostic verifies that the pinned fit changes the
source motion relative to convolution. These are numerical comparisons: the
pinned state-space fit has negative damping in active low-frequency joint
modes. Body surge, heave, pitch, PTO stroke, speed, force, and power are
checked at every saved time against the existing predeclared limits.
The external Applications/MATLAB-output path remains available for fresh
source runs. This fixture makes four already gated published cases durable in
ordinary Python CI; it does not add to the case badge count.

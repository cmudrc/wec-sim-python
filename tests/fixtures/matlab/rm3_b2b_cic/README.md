# RM3 body-to-body regularCIC published cases

`source.npz` contains the full 4,001-sample body and PTO CSV exports for
published `B2B_Case3` (without body-to-body hydrodynamics) and `B2B_Case4`
(with body-to-body hydrodynamics). They came from the successful
[pinned R2025b baseline run](https://github.com/cmudrc/wec-sim-python/actions/runs/38611588704)
on 2026-10-14, with WEC-Sim core commit
`0753b2e47f2457c078751dcfe5d251d1767b80ab` and Applications commit
`d53d4d4c9eda2581f04204f5d394a6ef84bb099e`. The archive was made by
loading the six corresponding source CSVs as float64 and calling
`numpy.savez_compressed` without rounding. Its SHA-256 is
`f27b5ab3397c8e7c112b398cfcc152f456c88270722e75ff827525b416005fc4`.

The portable gate runs each published 400 s input from the bundled RM3 HDF5,
with both the physical implicit added-mass scheme and the opt-in Simulink
delay scheme. It compares body surge, heave, pitch, PTO stroke, speed, force,
and power at every saved time with the existing predeclared numerical limits.
The external Applications/MATLAB-output path remains available for fresh
source runs. This fixture makes two already gated published cases durable in
ordinary Python CI; it does not add to the case badge count.

# RM3 fitted-radiation force histories

`source.npz` stores four unrounded float64 MATLAB body-output arrays for the
published 400 s `RM3_B2B_STATE_SPACE` Cases 5 and 6 (both bodies). Each array
has 4,001 rows and 37 columns. The archive was made directly from the four
corresponding CSV files in the successful [R2025b source run](https://github.com/cmudrc/wec-sim-python/actions/runs/37938077867),
with WEC-Sim core commit `0753b2e47f2457c078751dcfe5d251d1767b80ab`
and Applications commit `d53d4d4c9eda2581f04204f5d394a6ef84bb099e`.
Its SHA-256 is `a5d51a128d0b56651086a31f83010b16a3cb880363bec78587756b9e86509b55`.

The portable test uses the saved MATLAB velocities only as force-law inputs,
then compares independently replayed Python radiation forces with all six
saved MATLAB force channels on both bodies. The existing limit remains
`0.07%` of each channel peak, with a 3 N floor. This is a prescribed-motion
force check; the separate [RM3 B2B motion fixture](../rm3_b2b_cic/README.md)
checks the independent coupled trajectories. The pinned fitted radiation has
negative damping in active low-frequency joint modes, so this test checks
numerical agreement with that source fit, not its physical validity.

The MATLAB Applications input generates its `rm3.h5` through BEMIO and does
not track the HDF5. The replay passed with both the bundled HDF5 used in
ordinary Python CI and a locally regenerated HDF5 from the pinned BEMIO
input. The external Applications/CSV path remains available for fresh runs.

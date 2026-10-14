# Published WaveStar NMPC controller prefix

`controller_20s.npz` contains samples 0–400 (0–20 s at 0.05 s) from the
unchanged 225 s `WECCCOMP/WECCCOMP_Nonlinear_Model_Predictive` application.
The original [R2025b MATLAB run](https://github.com/cmudrc/wec-sim-python/actions/runs/37876181993)
used this repo at `11b4a264bf17a037cda9ad14b7c0d33f31ef54aa`, pinned WEC-Sim
core `0753b2e47f2457c078751dcfe5d251d1767b80ab`, and pinned Applications
`d53d4d4c9eda2581f04204f5d394a6ef84bb099e`. Its source script is
`tests/matlab/reference_model_baseline.m`, model `WECCCOMP_NMPC_SOURCE`.

The artifact `matlab-WECCCOMP_NMPC_SOURCE-baseline` supplies the MAT controller
trace and PTO CSV. Their SHA-256 hashes are, respectively,
`072ab1a3ee69ec97511d5a9225241845ffc73211213763e8a677644256a4612e`
and `4d09c2aaac0b8bf19c3520304abdcea1cf41e37a585f77acc7fa359e74a11b53`.
The fixture takes `time` and `command` from `cmd_ptoM`, `stroke` from
`motor_displacement`, `estimated_state` from `estimated_states`, and
`forecast` from `AR_excM_pred`. `axial_force` is PTO CSV column 33 (zero-based).
All arrays retain the source's float64 values. No source state or forecast is
fed into the Python observer, predictor, or controller during the test.

This fixture covers startup and five seconds of active NMPC. The published
full closed-loop plant motion remains open in `PARITY.md`.

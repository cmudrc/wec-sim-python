# WaveStar 225 s split-step source aggregate

`source.npz` is a lossless float64 subset of the full 225 s MATLAB R2025b
`WECCCOMP_NMPC_SPLIT_STEP_FULL_` export from
[run 39153812298](https://github.com/cmudrc/wec-sim-python/actions/runs/39153812298)
(`wavestar-nmpc-fine-active` artifact). The pinned source is WEC-Sim
`0753b2e47f2457c078751dcfe5d251d1767b80ab` and Applications
`d53d4d4c9eda2581f04204f5d394a6ef84bb099e`.

The saved arrays contain all 4,501 controller-clock samples from 0 to 225 s:
float pitch and pitch speed, PTO stroke and speed, axial force, requested
torque, and 500 sea components. They were copied without rounding from the
source body/PTO/component CSV files and `controller.mat`. The components
match the published `wavestar_nmpc_full_step/baseline.npz` array exactly.
`wavestar.h5` is the exact BEMIO file from the same artifact. Its 46 dataset
keys match `examples/data/wavestar.h5`, but eight fitted-radiation or IRF
datasets differ at floating-point precision. The largest absolute A-matrix
entry gap is `3.28e-12`; closed-loop feedback makes the exact source file
worth pinning.

SHA-256: `source.npz`
`be50ad71195b5ff2d20c97bb3270845354c5a703d3829bf4aa0e915196ba7866`;
`wavestar.h5`
`b2541ff01c319987d9157d1be9fbc9ec39a7660d0572d8e316c428b8aa3359d1`.

`test_wavestar_split_step_portable.py` independently advances the 1 ms
physical plant with the 50 ms controller and its one-sample source rate delay,
using only the pinned hydro data and incident sea phases. It checks both
fixture hashes before replay, and the macOS ARM64 and Linux x86 CI jobs run
the full portable test. It applies the
existing `<0.5%` limits to seven full-duration motion, force, command, and
PTO-work aggregates. The stricter pointwise autonomous gate remains failing:
current master differs by 0.007289 rad pitch against its 0.0005 rad limit
with this exact source HDF5.
This derived split-step gate does not qualify the 50 ms published case or
change the badge count.

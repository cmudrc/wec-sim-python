# WaveStar off-grid ode8 excitation

`source.npz` contains five un-ramped six-channel forces from the pinned
R2025b `irregExcF` function. The times are selected stage fractions of the
published 50 ms `ode8` clock, including samples between the 25 ms Python
excitation outputs. MATLAB used WEC-Sim `0753b2e47f2457c078751dcfe5d251d1767b80ab`
and the frequency coefficients interpolated from the archived
`../wavestar_fixed_input_step/wavestar.h5`; the 500 frequencies and phases
come from the published WaveStar NMPC input. SHA-256 of the fixture:
`52997c85e01f19b8439db37722f2b0101cc843da10b7367f23b8ac8ce53bb068`.

To regenerate, save `A`, `w`, `re`, `im`, `md`, `phase`, `dw`, `times`,
`direction`, and `spread` from the pinned components and HDF5 excitation
coefficients in a MAT file. Call
`tests/matlab/wavestar_ode8_stage_forcing_probe.m` with the pinned WEC-Sim
checkout, that MAT file, and an output CSV path. The source probe repeated
locally with exactly identical output. Python's direct component sum differs
from MATLAB by at most `4.98e-14` N or N m at these times; cubic
interpolation of the half-step force history differs by up to `3.63e-6`
N or N m. The portable test checks both raw and ramped stage forces.

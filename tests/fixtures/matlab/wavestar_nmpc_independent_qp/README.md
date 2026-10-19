# WaveStar NMPC QP on independently advanced full-duration inputs

`source.npz` contains all 4,201 controller samples from 15–225 s of a public
`run_wavestar_published(..., solver="wecsim_ode8")` Python solve. The solver
advanced its own WaveStar motion, observer states, and 40-step excitation
forecasts using only the pinned hydrodynamic HDF5 and 500 published wave
phases. The saved MATLAB plant motion, forces, controller states, and commands
were not supplied to that solve.

On 2026-10-18, local MATLAB R2025b Update 7 on mali ran the pinned
Applications commit `d53d4d4c9eda2581f04204f5d394a6ef84bb099e` `NMPC.m`
(SHA-256 `ca9a530630c46fefcd52a6f7791ecd523fa059027da91bd8921cdf70ae3157f4`)
with `quadprog` on those Python-generated observer states and forecasts.
The Python control path was from commit `bd3c124`; its WaveStar module is
unchanged from master at `c61a3f6`. The input MAT file had SHA-256
`de082b2daf69cec87ca5edd4e40a5f8e5a35cef83415276fb7edb1cbe827c9db`;
the direct MATLAB CSV had SHA-256
`ab5df0d9a11607de269c273b24ea910ea3761c6f0de98b076247cf8a3a285cc5`.

Maximum MATLAB/Python requested-torque difference is `0.000289571` N m,
RMS difference is `8.83404e-6` N m, and the 15–20 s maximum is
`1.02796e-5` N m. The independent run reaches 48 commands above 11.9 N m
in magnitude. These values are the archived L-BFGS-B controller outputs,
not the current solver's result. A local bounded least-squares replay on all
4,201 identical inputs reduces the maximum and RMS gaps to `0.000168290`
and `5.02630e-6` N m; its 15–20 s maximum is `5.20245e-6` N m. The
portable test now gates the recomputed commands at `0.00025` N m maximum,
`1e-5` N m RMS, and `1e-5` N m in the activation window.

On the independent full 225 s public source-`ode8` run, this QP change
reduces the local R2025b maximum pitch, speed, and axial PTO-force gaps
from `0.011410` rad, `0.056215` rad/s, and `4.725` N to `0.008283` rad,
`0.044096` rad/s, and `2.648` N. Relative signed PTO-work gap changes
from `0.1366%` to `0.05834%`. Both local runs pass their existing numerical
limits. Against the archived Linux source, the current local Python run still
misses the unchanged `0.012` rad pointwise pitch gate at `0.013908` rad.
These results do not establish portable physical-default parity or change
the published-case badge. A new local physical-default 225 s run with the
bounded least-squares controller still differs from the same R2025b source
by `0.041381` rad pitch, `0.443139` rad/s speed, `13.626` N PTO force, and
`5.9455%` signed PTO work, matching the previously documented physical gap.

To regenerate, run the public 225 s source-style solver with the published
HDF5 and wave phases, recording each `WaveStarNmpcObserver.step` and
`WaveStarNmpcPredictor.step` output. Save active samples as `states` (5×4,201),
`forecasts` (40×4,201), and `python_commands` (4,201×1), plus `time` and the
`A`, `Bu`, `Bw`, `R` matrices from
`../wavestar_nmpc_full_step/direct_controller_inputs.npz` in a MAT file. Then
run `tests/matlab/wavestar_nmpc_independent_qp_replay.m` with that MAT file and
the pinned Applications `WECCCOMP/WECCCOMP_Nonlinear_Model_Predictive`
directory. The direct output is time, MATLAB command, and Python command.

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
in magnitude. The portable test recomputes Python commands from all saved
independent inputs and gates their difference from the MATLAB `quadprog`
replay at `0.001` N m maximum, `0.00005` N m RMS, and `0.00005` N m in the
15–20 s activation window. These are controller-input component gates;
published autonomous closed-loop motion remains open.

To regenerate, run the public 225 s source-style solver with the published
HDF5 and wave phases, recording each `WaveStarNmpcObserver.step` and
`WaveStarNmpcPredictor.step` output. Save active samples as `states` (5×4,201),
`forecasts` (40×4,201), and `python_commands` (4,201×1), plus `time` and the
`A`, `Bu`, `Bw`, `R` matrices from
`../wavestar_nmpc_full_step/direct_controller_inputs.npz` in a MAT file. Then
run `tests/matlab/wavestar_nmpc_independent_qp_replay.m` with that MAT file and
the pinned Applications `WECCCOMP/WECCCOMP_Nonlinear_Model_Predictive`
directory. The direct output is time, MATLAB command, and Python command.

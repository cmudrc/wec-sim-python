# WaveStar NMPC QP with independent motion inputs

`source.npz` records 301 controller inputs from 15 to 30 s of the public
`run_wavestar_published(..., solver="wecsim_ode8")` path with an independent
`WaveStarNmpcPTO`. The Python run used only the pinned hydrodynamic HDF5 and
500 published wave phases. Its five observer states and 40-step excitation
forecasts, rather than MATLAB's saved states, were passed to the pinned
Applications `NMPC.m` under local MATLAB R2025b on mali on 2026-10-17. The
source function uses `quadprog` for each of the 301 active QPs.

The fixture stores the independent inputs, direct MATLAB commands, and
original Python commands. The maximum direct MATLAB/Python command gap was
`1.027964e-5` N m and RMS `6.053025e-7` N m. The portable test recomputes
Python commands from the saved independent inputs, then gates their difference
from the fresh MATLAB replay below `0.001` N m maximum and `0.00005` N m RMS.
It checks the optimizer on a trajectory off the saved MATLAB states; it does
not establish closed-loop motion parity.

To regenerate: run the public 30 s solver with the same published sea and
record each `WaveStarNmpcObserver.step` and `WaveStarNmpcPredictor.step`
output. Save the active samples (15–30 s) as `states` (5×301), `forecasts`
(40×301), and `python_commands` (301×1) in a MAT file. Copy `A`, `Bu`, `Bw`,
and `R` from `../wavestar_nmpc_full_step/direct_controller_inputs.npz` and
add `time` (301×1). In MATLAB, run
`wavestar_nmpc_independent_qp_replay(inputMat, applicationDir, outputCsv)`
using the pinned WaveStar NMPC application directory. The direct replay's
source is in `../../matlab/wavestar_nmpc_independent_qp_replay.m`.

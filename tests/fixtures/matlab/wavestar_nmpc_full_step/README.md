# Full-duration WaveStar NMPC timing diagnostic

`baseline.npz` retains the 225 s published 50 ms WaveStar NMPC trajectory from the first of two identical [R2025b source runs](https://github.com/cmudrc/wec-sim-python/actions/runs/38155596820). The workflow pinned WEC-Sim `0753b2e47f2457c078751dcfe5d251d1767b80ab` and WEC-Sim_Applications `d53d4d4c9eda2581f04204f5d394a6ef84bb099e`. SHA-256: `0205a5c8e7cb109bd6290a51ff6b42a5b69334c03fb2fc8ec445da3347aba6f6`.

The archive contains all 4,501 rows of float body motion and PTO output, the 500 wave components, 4,501 sea samples, and the controller command, power, energy, and their time grids. Arrays were loaded as float64 from the source CSV and MAT exports and compressed without rounding. The sea and first 200 samples of body and PTO output are exactly equal to the [independent fixed-input precontrol fixture](../wavestar_fixed_input_step/README.md). Its `wavestar.h5` has all 46 dataset values equal to this baseline's generated HDF5; the container SHA differs because HDF5 metadata differ.

`solver_25ms.npz` contains the full 225 s controller command, power, energy,
and their time grids from [R2025b run 38591012516](https://github.com/cmudrc/wec-sim-python/actions/runs/38591012516), plus the wave and component exports. The [follow-up R2025b run 38706528549](https://github.com/cmudrc/wec-sim-python/actions/runs/38706528549) reproduced those arrays exactly and added the float body and PTO histories at 50 and 25 ms respectively. The float is a 4,501-by-25 array; the PTO is a 9,001-by-43 array. All arrays were loaded as float64 without rounding and compressed. Updated SHA-256: `e847567a84a16ff709f98854fc801c45689ee204535c1609e7cdfa2589e0525e`. The source HDF5 has all 46 datasets equal to the earlier fixed-input fixture, including NaNs at the same positions; all 500 sea components and 4,501 wave samples are also equal. The first 9.95 s of body and PTO output reproduce the independent 25 ms precontrol fixture within `1e-9` in their units.

The refined run logs 9,001 controller samples at 25 ms instead of the published 4,501 at 50 ms. Between 30 and 225 s, the command changes on 3,899 of 3,900 intermediate half-step samples by more than `1e-6` N m. This means the attempted full-duration comparison changes controller behavior along with the solver step. The archived diagnostic proves that the full run cannot be interpreted as an isolated internal-plant-step test or a parity gate for the published 225 s case.

From the same pinned HDF5 and wave phases, an independent public Python 1 ms
plant run differs from the published 50 ms MATLAB trace by at most `0.04110`
rad pitch, `13.63` N PTO force, and `5.91%` PTO work. Relative to the 25 ms
MATLAB trace, the same Python run differs by `0.16693` rad, `70.66` N, and
`17.81%` work. The portable test checks this separation without claiming
published-case parity.

## Local MATLAB source run on mali

On 2026-10-17, R2025b on an Apple silicon Mac ran the unchanged published
225 s `ode8`/50 ms input after installing the Statistics and Machine Learning
and System Identification toolboxes. The source used the pinned WEC-Sim and
Applications commits above, copied the exact archived `wavestar.h5`, fixed
only the `.stl` filename case in a disposable application copy, and disabled
post-run plotting. `quadprog`, `dlqr`, `regress`, `ar`, and `forecast` all ran
from the installed MathWorks toolboxes. The source simulation took 337 s.

The local wave and 500 components match the hosted archive within `1e-12`.
The predeclared first-20-s body, PTO, force, and command gates passed. Through
10 s, float pitch differs from the hosted source by at most `4.1e-16` rad;
the source traces first differ by `1e-4` rad at 41.6 s. Over all 225 s, local
versus hosted source maxima are `0.007543` rad pitch, `0.05022` rad/s speed,
`1.507` mm PTO stroke, `2.858` N axial force, and `0.5884` N m command. The
closed-loop source is therefore sensitive to execution environment even with
identical published inputs and toolboxes.

An independent public Python 1 ms plant and 50 ms controller run against the
same waves differs from local MATLAB by at most `0.04138` rad pitch,
`0.4431` rad/s speed, `8.270` mm PTO stroke, `13.626` N axial force, and
`2.825` N m command. Its signed PTO work is `-159.666` J versus `-169.759` J
locally and `-169.702` J on the hosted run, a `5.95%` or `5.91%` relative
gap. The independent Python gap remains larger than the measured source
spread. This is diagnostic evidence; the full published case stays open.

`local_mali_audit.json` records every window and signal difference. Reproduce
it locally with `tests/matlab/wavestar_local_published_source.m` followed by
`tools/wavestar_local_source_audit.py --independent-python`. The MATLAB script
requires four arguments: the pinned WEC-Sim checkout, pinned Applications
checkout, this repository, and a fresh output directory. The Python script
takes that output directory as its positional argument.

## Reciprocal 225 s command and source-sensitivity checks

On 2026-10-18, local R2025b ran the source plant with the 4,201 active torque
requests from an autonomous Python `wecsim_ode8`/`WaveStarNmpcPTO` run. The
only intervention was replacing `y(1)` after the pinned `NMPC.m` computes
`y = Ubar`; MATLAB still advanced its own plant and evaluated the controller.
The source sea and all body/PTO samples before controller activation exactly
match an unchanged local source run. Source's logged active requests exactly
equal Python's saved requests. Full-duration source/Python maximum gaps are
`0.00013914` rad pitch, `0.00099627` rad/s speed, `27.36` µm PTO stroke, and
`0.002330` N active axial PTO force. Signed PTO work differs by `0.0520%`.
The [portable test](../../../test_wavestar_published_command_replay.py) uses
the same limits as the source-command plant gate, with the force check on the
active interval because source and Python have a small precontrol resistive
law difference. `python_command_source.npz` holds the source outputs and
Python requests; SHA-256 is
`79bf0453dd2e40987cf9a74433b89f75d748859387c81529a08a8fad1d154709`.
To regenerate the source run, export that fixture's active
`python_command` as a 4,201-by-1 `codexCommand` MAT array and pass its path as
the optional fifth argument to `wavestar_local_published_source`. The MATLAB
runner modifies only its disposable application copy.

A second local run changed one source NMPC request by `1e-4` N m on the 32nd
active controller call (16.55 s), after `y = Ubar`. The source wave and all
earlier motion are identical to the unchanged local run. This one-tick
perturbation grows to `0.0131005` rad maximum pitch difference at 120.55 s,
`0.0678542` rad/s maximum speed difference, and `1.19678` N m maximum later
request difference. `source_command_pulse.npz` stores the two complete pitch
and request histories; SHA-256 is
`50d7004c165bfd582c5a61687a9c144fdece57349897a0573dfb3c3412afa79a`.
The [source-sensitivity test](../../../test_wavestar_source_command_sensitivity.py)
checks the first changed tick and the subsequent amplification. This is a
MATLAB self-sensitivity diagnostic; it does not gate autonomous Python
motion. The published WaveStar NMPC case remains open.

A second local R2025b run changed only that same 16.55 s request by
`+0.00272438` N m, near the first independently advanced Python/source
request split. The original and changed source traces match exactly through
the command tick. Over 225 s, their maximum pitch and pitch-speed differences
are `0.00973867` rad and `0.0540725` rad/s. Peak absolute pitch changes by
`0.79430%`, while signed PTO work changes by `0.03518%`. The compact
[`source_command_pulse_2p7e.npz`](source_command_pulse_2p7e.npz) contains both
complete source histories; SHA-256 is
`c10290a6589c86fb8e4d49f71fd6739b9cd99f045062382708f4ac467f05c1de`.
The test above checks its first changed request and full-duration response.
The source inputs use the pinned WEC-Sim and Applications commits at the top
of this page. To regenerate, call `wavestar_local_published_source` twice
with fresh output directories: first with four arguments for the unchanged
case, then with an empty fifth argument and `.00272438` as the sixth.
This is a **modified published 50 ms source** control, not an independent
Python gate or the derived 1 ms plant/50 ms controller comparator. It does
not change the WaveStar case status or badge.

## Full-duration direct controller replay

`direct_controller_inputs.npz` retains all 4,201 active NMPC samples from the
same local R2025b source run, from 15 through 225 s. It contains the logged
five observer states, 40-step excitation forecast, and torque request, plus
the controller matrices used for an independent call to the pinned MATLAB
`NMPC.m` and `simulate_and_linearise.m`. The original source `controller.mat`
SHA-256 is `cf62ed39db7ddb9f0c0e9730bb60a9989a796a2c04d013bf4c800f4117a4c4b0`.
The pinned `NMPC.m` and `simulate_and_linearise.m` SHA-256 values are
`ca9a530630c46fefcd52a6f7791ecd523fa059027da91bd8921cdf70ae3157f4`
and `fa613a03240139303cd60edb43cd8838d703e3624702b06312f963def3862295`.
The compact fixture SHA-256 is
`3ff8f4ba8654ef66bf9058a9fe717cfa6c60868300d19e64a089fcb06781dbf0`.

The fresh direct MATLAB replay used the saved states and forecasts but
recomputed every QP with local `quadprog`; its 4,201 commands reproduce the
logged source within `4.39e-11` N m. The exported direct-command CSV SHA-256
is `e3e6bc5ddb902c3dd5c616b9a2bbb002d09af6dc0cfb398d1c9bf81f728d0160`.
To repeat it, save the fixture's `A`, `Bu`, `Bw`, `R`, transposed
`estimated_state`, `forecast`, and column `source_command` as `A`, `Bu`,
`Bw`, `R`, `states`, `forecasts`, and `source` in a MAT file; then call
`wavestar_nmpc_direct_replay(applicationDir, inputMat, outputCsv)` from
`tests/matlab` with the pinned application function directory.

Given these saved source inputs, independent Python `WaveStarNmpcController`
commands differ from the fresh direct MATLAB replay by at most `0.0007652`
N m and `0.0000203` N m RMS over all 225 s. The CI limits are `0.001` N m
maximum and `0.00005` N m RMS, plus `5e-6` N m maximum through 20 s. This
is a direct-input controller gate. It does not advance the Python observer,
predictor, or plant, and does not close published WaveStar motion parity.

![Full-duration direct controller torque and gap](../../../../docs/figures/wavestar-nmpc-direct-controller-gap.png)

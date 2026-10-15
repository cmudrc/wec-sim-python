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

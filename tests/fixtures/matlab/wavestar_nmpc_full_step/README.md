# Published WaveStar NMPC baseline for the internal-step audit

`baseline.npz` retains the 225 s published 50 ms WaveStar NMPC trajectory from the first of two identical [R2025b source runs](https://github.com/cmudrc/wec-sim-python/actions/runs/38155596820). The workflow pinned WEC-Sim `0753b2e47f2457c078751dcfe5d251d1767b80ab` and WEC-Sim_Applications `d53d4d4c9eda2581f04204f5d394a6ef84bb099e`. SHA-256: `0205a5c8e7cb109bd6290a51ff6b42a5b69334c03fb2fc8ec445da3347aba6f6`.

The archive contains all 4,501 rows of float body motion and PTO output, the 500 wave components, 4,501 sea samples, and the controller command, power, energy, and their time grids. Arrays were loaded as float64 from the source CSV and MAT exports and compressed without rounding. The sea and first 200 samples of body and PTO output are exactly equal to the [independent fixed-input precontrol fixture](../wavestar_fixed_input_step/README.md). Its `wavestar.h5` has all 46 dataset values equal to this baseline's generated HDF5; the container SHA differs because HDF5 metadata differ.

The full-duration audit changes only the pinned Simulink `ode8` internal step to 25 ms. Its MATLAB controller, WEC-Sim, wave, and body-output grids must remain at 50 ms. The archive supplies the unchanged published-step comparator; it is not a new parity gate for the published 225 s controlled case.

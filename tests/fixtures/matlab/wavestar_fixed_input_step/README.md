# WaveStar fixed-input internal-step fixture

These arrays come from the passing [R2025b source run](https://github.com/cmudrc/wec-sim-python/actions/runs/38583530934) of PR #357. The workflow checks out WEC-Sim `0753b2e47f2457c078751dcfe5d251d1767b80ab` and WEC-Sim_Applications `d53d4d4c9eda2581f04204f5d394a6ef84bb099e`. It holds the published 50 ms wave, WEC-Sim, controller, and body-output settings fixed, and runs the precontrol interval through 9.95 s at Simulink ode8 internal steps of 50, 25, and 12.5 ms. The PTO logger follows the internal step.

`source.npz` contains the body-1 and PTO-1 records, 500 JONSWAP components, wave history, and step settings from all three runs, compressed without reducing precision. `wavestar.h5` is the exact BEMIO output from that source run. SHA-256 values:

| File | SHA-256 |
| --- | --- |
| `source.npz` | `f6c0b2b189d7ae4603ed24eb24d6c312f9e35bd90310d8e40fa0a2a7c476ade6` |
| `wavestar.h5` | `22702b94f052dfe6bb93d257006654d4e51246020a53f25b432d591bf8139d1a` |

`tests/test_wavestar_nmpc_internal_step.py` uses this fixture in ordinary Python CI. Setting `WEC_SIM_WAVESTAR_INTERNAL_STEP_DIR` and `WEC_SIM_WAVESTAR_H5` makes the same test read fresh MATLAB exports instead. The gate checks source-input identity, zero PTO actuation, output grids, and the independent Python pitch comparisons. It does not claim parity for the published 225 s controlled case.

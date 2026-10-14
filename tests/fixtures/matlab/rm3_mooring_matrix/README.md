# Published RM3 MooringMatrix, 400 s

`source.npz` retains all 40,001 samples of the wave, active surge/heave/pitch
positions and speeds of both bodies, PTO stroke/speed/force/power, and mooring
surge position/force from the [passing R2025b source run](https://github.com/cmudrc/wec-sim-python/actions/runs/37702732873).
It is a lossless float64 selection of that run's CSV columns, without rounding.
SHA-256: `cd0bae84c6c70f95fcbf461b8cb800391e233742ca5531bd4b60a28d6e6e65f9`.

`etaData.mat` is the published input from WEC-Sim_Applications commit
`d53d4d4c9eda2581f04204f5d394a6ef84bb099e` at
`Mooring/MooringMatrix/etaData.mat` (SHA-256
`73c5c931d74f5f956c035262e1b4c44406603d688d59300e2363fbda31110979`).
`rm3.h5` is the R2025b BEMIO-generated HDF5 from the same source run
(SHA-256 `af41df46b066e89c70cf931b9aeaf33be07f0a0aaf7f111dacfd5679ef9807a6`).

The ordinary Python CI test advances the public `WEC.run` floating joint
independently from these pinned inputs for the full 400 s and applies the
previously declared body/PTO/mooring limits. It also checks signed PTO work
within 1% of the source. No saved source motion or forces drive the run.

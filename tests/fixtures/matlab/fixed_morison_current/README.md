# Fixed monopile current source

`source.npz` stores the full 20 s, 10 ms output from the three R2025b
`MORISON_FIXED_CURRENT` cases in [run 37820035413](https://github.com/cmudrc/wec-sim-python/actions/runs/37820035413),
artifact `matlab-MORISON_FIXED_CURRENT-baseline`. The workflow used pinned
WEC-Sim core `0753b2e47f2457c078751dcfe5d251d1767b80ab` and Applications
`d53d4d4c9eda2581f04204f5d394a6ef84bb099e`; its repository commit was
`7d84c79269b51ed13002962045041f0b0242787e`.

The shared 500-component phase table, 2,001-sample wave, and stationary tower
trace are stored once. Each current profile has its own complete six-component
source force and monopile trace. Arrays were loaded from the artifact CSVs
with `numpy.loadtxt` and saved losslessly with `numpy.savez_compressed`.
The archive SHA-256 is
`203e0e7e702655b18c98b47bbc32ea2cab4f447c4859a0f95a34b4f75d972c07`.
The three source force CSV hashes, in uniform, power, and linear order, are
`e4d7360ae958065e6e966241ac26a0a222a64e08c48434418df3b79a18746b50`,
`030db01942105c159b4ae2f0aaeafa12923311ce3e69b7ddacdb799ef641cb4a`,
and `2f006f212112916e4889b11c0667f0a17feae3188c6993cf986fb323edb68e55`.

The ordinary Python test independently advances the published PM sea and
fixed Morison element for every sample. It uses saved phases as the stochastic
input, and checks the wave within `1e-10` m and each force or moment within
`1e-3` N or N m. MATLAB body force and stationary-body histories are also
checked. No MATLAB force or body trajectory is advanced as a Python input.

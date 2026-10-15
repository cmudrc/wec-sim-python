# Published OSWEC Full_Directional_Waves, 400 s

`fullDirSpectrum.mat` is the published spectrum from WEC-Sim_Applications
commit `d53d4d4c9eda2581f04204f5d394a6ef84bb099e` (SHA-256
`c0abee4fe5da352412d64342ea276832bfb19595774a0f1f2106f5f109dfaee9`).

`source.npz` selects unrounded float64 columns from the [passing R2025b
source run](https://github.com/cmudrc/wec-sim-python/actions/runs/38663378748):
all 8,001 wave samples, the 47-by-180 phase bank, both bodies' six-position,
six-velocity, and six-excitation histories, plus PTO torque and power. SHA-256:
`52a83e9d86b18b16ea876d6d8ea88cd6c04db80b189ba1717d3b661fac6e121b`.

The portable gate reuses `../oswec_multi_wave/oswec.h5` (SHA-256
`b9ccf3a32e0893c23b1645c15a2baef4a0c6f678b90d177d5797e3c65de4b7d3`).
An array-by-array audit against the source run's BEMIO HDF5 found identical
frequency-domain coefficients; only four time-domain impulse-response arrays
differ, by at most `4.1e-11` in stored units. This avoids adding a second large
OSWEC HDF5. The public `WEC.run` gate advances the full 400 s independently;
neither source motion nor source excitation drives the result.

The source force block omits heading-bin width. The gate requests that
source-compatible option explicitly and checks the physically integrated force
separately. It does not endorse the omitted-width force as a physical default.

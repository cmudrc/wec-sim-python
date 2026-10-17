# RM3 small-wave source record

Generated locally with MATLAB R2025b and pinned WEC-Sim source
`0753b2e47f2457c078751dcfe5d251d1767b80ab`. The input is
[`examples/data/rm3.h5`](../../../../examples/data/rm3.h5)
(SHA-256 `3bdfdcfa3dcd2db1851ab2026eb6b4fdc565ad86be8503e09b42231feee429fb`).

[`rm3_small_wave_scaling_reference.m`](../../../matlab/rm3_small_wave_scaling_reference.m)
generates two 20 s PM incident seas at heights 0.01 and 0.02 m with
32 equal-energy components, substream seed 4, a 5 s ramp, and 0.1 s samples.
Set `WEC_SIM_SOURCE_DIR` to the pinned source directory and
`WEC_SIM_MATLAB_RM3_SCALING_DIR` to an output directory, then run the script
in R2025b. The CSVs retain MATLAB's spectral bins, phases, and elevation;
the ordinary Python test independently rebuilds these and advances RM3
motion and PTO work at both heights.

| File | SHA-256 |
| --- | --- |
| `phase_seed4.csv` | `b3ccba29969002a83f580781c2712ae8f24a8835cb32f8958076214cb034ffe2` |
| `bins_001.csv` | `c4fc199164a072eb25a63e4980e09fbebde0259e0d4550b84909d33abb8859c3` |
| `bins_002.csv` | `141843baf59278034b3cabc797da9c4ddcc0e9d879f2c623c69304e5c6465ba3` |
| `wave_001.csv` | `a7b97e0ded990d5c922f6ef7b93f77aa685953f109f46d27239834633429704e` |
| `wave_002.csv` | `c1404d1505bb7556313150ab1397e679a32c01336c80a09e04f6ffce254a655a` |

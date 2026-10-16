# Floating OWC published wave and body-force fixture

The [pinned R2025b published `OWC/FloatingOWC` run](https://github.com/cmudrc/wec-sim-python/actions/runs/37895804527)
produced `matlab-FLOATING_OWC_SOURCE-baseline`. This folder contains its
`applications/OWC/FloatingOWC/hydroData/floatingOWC.h5` unchanged and a
compressed subset of three output CSV files. The HDF5 SHA-256 is
`c9092b10cf35bc438805f719da5434478a2211ca9261e9a33d562dbf4d4e1a11`.

`source.npz` stores the common 50,001-sample time grid, wave elevation, and
the six excitation-force channels from columns 20–25 (one-based) of each
published body CSV. During extraction, both body time columns matched the wave
time column exactly. Original CSV SHA-256 values:

| Source file | SHA-256 |
| --- | --- |
| `FLOATING_OWC_SOURCE_wave.csv` | `576845e62a1f613e4329838317f4e760f458cafe5d91fb93ec8a924200a3f700` |
| `FLOATING_OWC_SOURCE_FloatingOWC_body1.csv` | `224353c7b13f6cc0985a327cff8ab217f503d5c447fcf61d6928936fb069b699` |
| `FLOATING_OWC_SOURCE_FloatingOWC_body2.csv` | `1a4730dfc74e6e478e0ac66c28d24a49cbf6a60c6cb3aa7935643b2ec43b711f` |

The [portable test](../../../test_floating_owc_wave_force_portable.py)
recomputes the regular wave and both bodies' excitation from the exact HDF5
and published input settings. It does not use the saved body states, mooring
forces, or chamber states. The full coupled-motion gate remains in the MATLAB
workflow; this fixture makes its wave and hydrodynamic-force inputs checkable
on every Python CI run without native MoorDyn.

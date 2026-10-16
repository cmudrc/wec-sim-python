# Floating OWC published wave and body-force fixture

The [pinned R2025b published `OWC/FloatingOWC` run](https://github.com/cmudrc/wec-sim-python/actions/runs/37904325081)
produced `matlab-FLOATING_OWC_SOURCE-baseline`. This folder contains its
`applications/OWC/FloatingOWC/hydroData/floatingOWC.h5` unchanged and a
compressed subset of eleven output CSV files. The HDF5 SHA-256 is
`f8d225dcb039850616f444e25bf0c3a93212ffce82508471dcb517bbbde49bcc`.

`source.npz` stores the common 50,001-sample time grid, wave elevation, and
the six excitation-force channels from columns 20–25 (one-based) of each
published body CSV. It also stores each body's logged rigid mass and inertia,
6×6 added mass, radiation damping, and hydrostatic matrices. During extraction,
both body time columns matched the wave time column exactly. Original CSV
SHA-256 values:

| Source file | SHA-256 |
| --- | --- |
| `FLOATING_OWC_SOURCE_wave.csv` | `576845e62a1f613e4329838317f4e760f458cafe5d91fb93ec8a924200a3f700` |
| `FLOATING_OWC_SOURCE_FloatingOWC_body1.csv` | `224353c7b13f6cc0985a327cff8ab217f503d5c447fcf61d6928936fb069b699` |
| `FLOATING_OWC_SOURCE_FloatingOWC_body2.csv` | `1a4730dfc74e6e478e0ac66c28d24a49cbf6a60c6cb3aa7935643b2ec43b711f` |
| `FLOATING_OWC_SOURCE_body1_mass.csv` | `e40dca224f0d3a8dda1944569e7a68ab94b0231ba6cb409abd3f7bc4d6e69207` |
| `FLOATING_OWC_SOURCE_body1_added_mass.csv` | `1ec279fd0221d32e1a68cef166a1a2bed6d0c1043d907c77cbe2ecc3daf4a6c0` |
| `FLOATING_OWC_SOURCE_body1_radiation_damping.csv` | `ee794dce1d701fdc314dfdd6ea015af085ca016c7153307da7dd7fe0abc74676` |
| `FLOATING_OWC_SOURCE_body1_hydrostatic.csv` | `2d334516393cacfc1df11e094574261e7da39fca576ba44d5edd3f64d2b164c3` |
| `FLOATING_OWC_SOURCE_body2_mass.csv` | `01f2edb1a242d211f60f626051a38582072cc409a8daec549f1c0bf5308a0ec4` |
| `FLOATING_OWC_SOURCE_body2_added_mass.csv` | `1e17243e547d45e3df5dad9d1fb1e23306741a76bdd32b5c8ac169c24c5983e3` |
| `FLOATING_OWC_SOURCE_body2_radiation_damping.csv` | `ea2796497a6daf9c1786ba9a582d47c8d786f51cd75a7199975c81f096b4066b` |
| `FLOATING_OWC_SOURCE_body2_hydrostatic.csv` | `73ec268ecc9e7858c2d7b0dfb1e75918ac677d69aa9efbc8b16102fd8202dea4` |

The [portable test](../../../test_floating_owc_wave_force_portable.py)
recomputes the regular wave, both bodies' excitation, and hydrodynamic matrices
from the exact HDF5 and published input settings. It does not use the saved body
states, mooring
forces, or chamber states. The full coupled-motion gate remains in the MATLAB
workflow; this fixture makes its wave and hydrodynamic-force inputs checkable
on every Python CI run without native MoorDyn.

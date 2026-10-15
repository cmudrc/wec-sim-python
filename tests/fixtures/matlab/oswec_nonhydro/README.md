# Published OSWEC fixed-base reference

`source.npz` contains the full 4,001-sample flap, fixed base, and PTO output
from the pinned R2025b [`OSWEC_Nonhydro` source run](https://github.com/cmudrc/wec-sim-python/actions/runs/37795143695).
That run uses WEC-Sim core `0753b2e47f2457c078751dcfe5d251d1767b80ab`
and Applications `d53d4d4c9eda2581f04204f5d394a6ef84bb099e`. Each NPZ
array retains its source CSV columns: time, six positions, six velocities,
six forces, and six excitation or power components.

The source CSV SHA-256 values are:

- Flap: `acc0bbff9ccd4c803529193a7d9b47cd175c688a1da71a8c536aa33a5ff9fd24`
- Base: `5a4d4d24b77a570cc76f487a0f0288ec02a4148c3fae1c1dbb761925e0083c43`
- PTO: `456f7371b7a3ecc3e95ec0dc30f9b5c722650bde1897f02e8b817f28e33c3894`

The NPZ SHA-256 is
`990ca6819499d2769f44bd3a22cfaca7e4b9c32093ab960dc3d4b18eba55dfa4`.
The test reuses the 11 MB BEMIO OSWEC HDF5 already archived for the multiple
wave case (`b9ccf3a32e0893c23b1645c15a2baef4a0c6f678b90d177d5797e3c65de4b7d3`).
All 68 datasets have the same shapes and names as this run's generated HDF5;
the largest numeric difference is `4.01e-11` in the radiation IRF, with no
other difference above `1e-9`.

`tests/test_oswec_nonhydro_parity.py` still accepts live source files when
provided, and otherwise independently advances the full published 400 s case
against this archived trace in ordinary Python CI. Its existing position,
velocity, excitation, and PTO limits are unchanged. This strengthens a
previously counted case; it does not add a badge gate or compute base reactions.

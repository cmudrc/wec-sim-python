# Published RM3 radiation-options fixture

The successful [R2025b reference sweep](https://github.com/cmudrc/wec-sim-python/actions/runs/38165906178)
uploaded `matlab-RM3_Radiation_Options-baseline` at run commit
`e6330e5eb34529f115f4d0be3f894b563f8984a3`. Its workflow pins WEC-Sim
core `0753b2e47f2457c078751dcfe5d251d1767b80ab` and Applications
`d53d4d4c9eda2581f04204f5d394a6ef84bb099e`. The published application
ran all four radiation settings for 500 s with a 2.5 m, 12 s regular wave.

`rm3.h5` is the exact BEMIO-generated HDF5 from that artifact (SHA-256
`75e87c835ceb06984e425030bfef2e4c29fc0c6c63a73a1b140e62cb4828d6a9`).
It differs from `examples/data/rm3.h5`, so the bundled file cannot replace
this published-case input.

`source.npz` contains all twelve source CSV arrays as float64 without
rounding: two 5,001-sample body histories and one PTO history for each of
constant, convolution, FIR, and fitted state-space radiation. Every array
was checked for exact numerical equality with its downloaded CSV after
`numpy.loadtxt`; the archive SHA-256 is
`426c9e1a233fc6f52520f29d4fcb452a06edc1ed3d0423f4911521e2b76ac8e6`.

`tests/test_rm3_radiation_options_parity.py` now runs its existing independent
full-duration motion, PTO, and radiation-force gates in ordinary Python CI
from this fixture. Setting both `WEC_SIM_APPLICATIONS_DIR` and
`WEC_SIM_MATLAB_MODEL_OUTPUT_DIR` runs the same tests against fresh MATLAB
outputs. This preserves gates for an already counted RM3 case; it does not
add a case badge or establish passivity of the source's fitted radiation.

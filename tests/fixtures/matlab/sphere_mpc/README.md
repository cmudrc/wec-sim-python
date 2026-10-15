# Published Sphere Controls/MPC, 400 s

The [passing R2025b baseline run](https://github.com/cmudrc/wec-sim-python/actions/runs/37743807105)
uses WEC-Sim_Applications commit `d53d4d4c9eda2581f04204f5d394a6ef84bb099e`
and core WEC-Sim commit `0753b2e`. `sphere.h5` is its BEMIO-generated HDF5;
`coeff.mat` is the published fourth-order radiation fit; `setup.mat` contains
the exported MATLAB prediction and QP matrices.

`source.npz` retains all 40,001 exported samples of wave elevation, Sphere
heave position and speed, heave excitation, PTO force and power, and internal
MPC plant speed and position. It also retains the 500 JONSWAP components,
all 801 force-rate commands, and every 31-value source excitation forecast.
The arrays are selected from the exported MATLAB CSV without interpolation or
additional rounding. The source motion, excitation, and command histories
are comparison targets only; the Python closed loop advances from the pinned
HDF5, coefficient file, and wave phases.

SHA-256 digests are asserted in the portable test for all four files. The
fresh hosted MATLAB gate separately checks the full source export and
component-level behavior. Ordinary Python CI runs the independent full 400 s
closed loop against this pinned selection.

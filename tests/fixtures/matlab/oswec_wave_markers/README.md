# Published OSWEC Wave Markers hydrodynamic input

The pinned MATLAB Applications `Wave_Markers/OSWEC` workflow generated
`oswec.h5` by running `_Common_Input_Files/OSWEC/hydroData/bemio.m` against its
checked-in `oswec.out` with WEC-Sim core
`0753b2e47f2457c078751dcfe5d251d1767b80ab` and Applications
`d53d4d4c9eda2581f04204f5d394a6ef84bb099e`. The archived source file
has SHA-256
`79d9ab232958f3b3f9ec1bcdcd9a56c9fe8eb54f22e9de2714403860b8a51b32`.

This fixture is the same source HDF5 repacked with gzip and shuffle filters to
keep the portable gate small. All 68 datasets retain their original values,
shapes, and attributes, including the 36 wave headings. Its SHA-256 is
`f99423db63bda6b744bec8ca88e06baedd3475b76c8a2282679c46f18907b374`.
The published 400 s trajectory is in
[`../oswec_wave_markers.npz`](../oswec_wave_markers.npz); the portable gate is
[`../../../test_oswec_wave_markers_source_ode4.py`](../../../test_oswec_wave_markers_source_ode4.py).

The repo's `examples/data/oswec.h5` is a different 10-heading file. It cannot
serve as the exact hydrodynamic input for this published application.

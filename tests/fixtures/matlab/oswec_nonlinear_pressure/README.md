# OSWEC nonlinear visualization pressure samples

`samples.npz` contains 46 unrounded float64 frame samples from the untouched
published 120 s R2025b
[MATLAB run](https://github.com/cmudrc/wec-sim-python/actions/runs/38259283403)
(`matlab-oswec-nonlinear-viz-pressure`, run commit
`1083a1f29725785cbba0cb67f8e389a9cade16fc`). The workflow pins WEC-Sim
core `0753b2e47f2457c078751dcfe5d251d1767b80ab` and Applications
`d53d4d4c9eda2581f04204f5d394a6ef84bb099e`. The original `source.mat`
SHA-256 is `f5c2944777f35741a48034e64fb4d3d7ec165466f6f2fc76f7be1f38383ad466`.

The selected zero-based frame indices are every 30 samples from 0 through
1200, plus 1, 10, 399, 400, and 401. This spans startup, the 40 s ramp
boundary, and the full duration. Load `source.mat` with `scipy.io.loadmat`,
select these rows from `time`, `pose`, `wave`, `hydrostatic`, `nonlinearWave`,
`linearWave`, and `forceRestoring`, and copy `cg`, `centers`, `normals`, `areas`,
`waterDepth`, `waveNumber`, `rho`, and `gravity` without interpolation or
rounding. Load the pinned Applications
`_Common_Input_Files/OSWEC/geometry/flap.stl` with
`trimesh.load_mesh(..., process=False)` and save its `vertices` and `faces`
arrays in the same archive with `numpy.savez_compressed`. The archive SHA-256 is
`a6b808a18fe6b2e7cfe6a6e4abdf45ff0c7bdb51f29d16752414e5a6c32a0005`.

The source STL SHA-256 is
`ed171c2bf7de8162b048ef1f1498bc1ee5200315c3c1d91b6b319fd5f1d15769`.
Ordinary Python CI compares all 1,042 facet pressures in each selected frame
and their integrated hydrostatic restoring forces. The dedicated MATLAB
workflow still checks every saved frame. Source poses drive both pressure
checks; this fixture does not establish independent published-step motion
parity or change a case badge.

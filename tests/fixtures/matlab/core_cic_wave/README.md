# Core CIC wave preprocessing

These compact CSVs were generated locally with MATLAB R2025b Update 7 from
pinned WEC-Sim source commit `0753b2e47f2457c078751dcfe5d251d1767b80ab`
using [`tests/matlab/core_cic_wave_reference.m`](../../../matlab/core_cic_wave_reference.m).
The input uses finite 25 m depth, BEM frequencies 0.4–2 rad/s, and a 2 s
time grid. `noWaveCIC` leaves period unspecified so the source selects the
minimum BEM frequency. `regularCIC` uses a 30° heading and three off-origin
markers.

Regenerate from the repository root after adding the pinned MATLAB source
tree to the path:

```matlab
addpath(genpath('/path/to/pinned/WEC-Sim/source'))
run('tests/matlab/core_cic_wave_reference.m')
```

The paired Python test checks all samples, the selected frequency and period,
finite-depth wave number, amplitude, marker elevations, and regular-wave
power. This checks wave preprocessing, not RM3 body-motion parity; published
case badges do not change.

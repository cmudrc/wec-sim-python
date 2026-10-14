# Directional irregular-wave gauge probe

These CSVs were generated locally on 2026-10-14 with MATLAB R2025b Update 7
and pinned WEC-Sim source commit `0753b2e47f2457c078751dcfe5d251d1767b80ab`.
They retain 16 traditional PM frequency bins, the 16-by-2 phase matrix, and
nine 50 ms samples of origin and three marker elevations. This is a wave-input
component comparison, not a full published-case motion gate.

After adding the pinned WEC-Sim `source` directory to the MATLAB path:

```matlab
w = waveClass('irregular');
w.height = 2;
w.period = 8;
w.direction = [0,45];
w.spread = [0.6,0.4];
w.spectrumType = 'PM';
w.bem.option = 'Traditional';
w.bem.count = 16;
w.phaseSeed = 3;
w.marker.location = [3,4; -2,1; 0,0];
time = (0:0.05:0.4)';
w.setup([0.2,2], 15, 0.2, 0.05, 8, time, 9.81, 1025);
writematrix([w.waveAmpTime,w.waveAmpTimeViz(:,2:4)], 'elevation.csv');
writematrix([w.omega,w.dOmega,w.wavenumber,w.spectrum,w.amplitude], 'components.csv');
writematrix(w.phase, 'phase.csv');
```

Column order in `elevation.csv` is time, origin, marker `[3,4]`, marker
`[-2,1]`, and marker `[0,0]`. Maximum Python/source differences were
`8.33e-16` m for elevation, `7.77e-16` across component values, and
`5.33e-15` rad for phase.

Both languages reproduce the pinned source's finite-depth wavenumber
iteration. At the first `0.2` rad/s bin, its dispersion residual is `2.76%`;
that bin's PM amplitude is only `2.35e-126` m²s/rad. The energy-weighted
relative residual over all bins is `3.83e-14`. The test records this source
limitation and checks dispersion on bins carrying appreciable energy.

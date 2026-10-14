# Finite-depth regular-wave gauge probe

These CSVs were generated locally on 2026-10-14 with MATLAB R2025b Update 7
and pinned WEC-Sim source commit `0753b2e47f2457c078751dcfe5d251d1767b80ab`.
They contain nine 50 ms samples of origin and three marker elevations, plus
finite-depth wavenumber and wave power. The Python test also checks the
dispersion relation and directional phase analytically. This is a wave-input
component check, not a published-case motion gate.

After adding the pinned WEC-Sim `source` directory to the MATLAB path:

```matlab
w = waveClass('regular');
w.height = 2;
w.period = 8;
w.direction = 30;
w.marker.location = [3,4; -2,1; 0,0];
time = (0:0.05:0.4)';
w.setup([0.1,2], 15, 0.2, 0.05, 8, time, 9.81, 1025);
writematrix([w.waveAmpTime,w.waveAmpTimeViz(:,2:4)], 'elevation.csv');
writematrix([w.wavenumber,w.power], 'constants.csv');
```

Column order in `elevation.csv` is time, origin, marker `[3,4]`, marker
`[-2,1]`, and marker `[0,0]`. The observed maximum Python/source elevation
gap was `4.44e-16` m, wavenumber gap `2.78e-17` rad/m, and power gap
`1.46e-11` W/m.

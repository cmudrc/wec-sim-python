# Seeded two-column imported spectrum

`spectrum.mat` contains seven frequency/density rows and no phase column. The
first and last rows lie outside the bundled Sphere HDF5 BEM range. Pinned
WEC-Sim `0753b2e47f2457c078751dcfe5d251d1767b80ab` generated
`bins_phase.csv` and `elevation.csv` with MATLAB R2025b Update 7,
`waves.phaseSeed = 1`, a 0.1 s grid through 20 s, and a 1 s ramp. The
reproducible source command is
`tests/matlab/imported_spectrum_two_column_reference.m`.

SHA-256: `spectrum.mat` `1d79fbb0a4855fddbf9906b7c9e16307a6359c191bee9b512da6964aecba9ea2`;
`bins_phase.csv` `b556a831b75c153739c991b3bd9ef272887bc96a314acde19fda42f9d3b97533`;
`elevation.csv` `328147d68a2d8ad5dcde5afff019cf9170780608f833484187840e03ed3ac5b9`.

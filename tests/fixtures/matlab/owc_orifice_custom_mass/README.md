# OWC orifice custom-mass source check

These 11-sample CSVs came from local MATLAB R2025b Update 7 using pinned
WEC-Sim core `0753b2e47f2457c078751dcfe5d251d1767b80ab` and Applications
`d53d4d4c9eda2581f04204f5d394a6ef84bb099e`. The starting model was
`OWC/OrificeModel/OWC_GBM.slx`, its published input, geometry, orifice data,
and the checked-in `../owc_orifice.h5` (SHA-256
`fe702ecf99f10d4cd0bba5d858b588dd37709cd45a5c3727356849e8dc40d676`).

For this bounded diagnostic, set `simu.endTime = 0.05` s and
`waves.height = 1e-8` m. The heavy run sets `body(1).mass = 689.049` kg,
100 kg above displaced-water mass; the control keeps `mass = 'equilibrium'`.
All other physical settings, including `ode23t`, the 0.005 s output step,
orifice, damping, drag, and 15 s radiation memory, remain the published
values. The local installation lacks the Control System Toolbox. The input's
unused `tf` filter declarations were removed for both runs; `OWC_GBM.slx`
contains no reference to those filter variables. The tiny wave leaves the
equilibrium control below `3e-16` m heave by 0.05 s.

After `wecSim`, export each run from its MATLAB workspace:

```matlab
values = [output.bodies(1).time(:), output.bodies(1).position(:,3), ...
          output.bodies(1).velocity(:,3)];
writematrix(values, 'heave.csv');
```

The CSV columns are time, heave, and heave speed. The heavy `heave.csv` SHA-256
is `6ea928ba64ff66d7acb131606cb94dff4b4845f0e405674b7e4fbcea83909745`;
the control `equilibrium_heave.csv` SHA-256 is
`67f14627c37e4670a30fa8b9c75074f16affdcfb5c5705b7a85804f7a736fdfd`.
An independent public Python run with the explicit source-convention orifice
route differs by at most 4.48 µm heave and 0.522 mm/s heave speed. This is a
custom-mass diagnostic, not a new full-duration published-case gate.

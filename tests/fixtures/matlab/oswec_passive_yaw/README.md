# Published OSWEC passive yaw off and on

`off.npz` and `on.npz` retain the active 60,001-sample motion, PTO, and wave
histories from the published 600 s `Passive_Yaw/PassiveYawOFF` and
`Passive_Yaw/PassiveYawON` cases. The ON fixture also retains all six source
excitation channels. They were packed without rounding from the successful
[R2025b baseline run](https://github.com/cmudrc/wec-sim-python/actions/runs/38677803196)
(run commit `48e11f79682346404acb3762fc80d10ddcf8a34c`). Their SHA-256 values
are `6a567ae9d52b54f6abf1505629049f64edfa13d7b8786366ecccfc07d77f7d7e`
and `d51699a5ea3b1cf5c8894aacfda9665618617d9b770734317b526dcd33590295`.
The original ON and OFF body-1 CSV SHA-256 values are
`b792e9cd28b20b5901dd90422cac432a190d287778d621d08728c4044725689b`
and `276e5f7a4e5035ea612ec0455a36616b946bfb9a42de21e06e8c7c9c5d82536a`.
To reconstruct either fixture from that run's `OSWEC_PASSIVE_YAW_*` CSVs,
take time, yaw, yaw speed, and (ON only) six excitation components from
body 1 columns 0, 6, 12, and 19:25. Take PTO angle, speed, torque, and
power from PTO columns 5, 11, 17, and 23; wave elevation is wave column 1.
The stationary flap coordinates are body 1 columns 1:6 and 7:12; the base
position and velocity are body 2 columns 1:7 and 7:13.

The four inactive-coordinate histories are represented by their first sample
and measured maximum departure from that sample. The largest source
departure is `6.04e-16`; this keeps the portable arrays compact while the
test still checks every Python sample against the source pose and the
recorded source drift.

The source HDF5 SHA-256 is
`9e5570674735bc62a1272077a08110fda91262e38769c900cc33eb4dfa0f5970`.
The bundled `oswec_multi_wave/oswec.h5` has the same 68 dataset paths;
all but four impulse-response arrays are exactly equal, and the largest
absolute difference in those four is `4.01e-11`. Using the bundled file,
the full-duration local OFF maximum yaw gap is `3.65e-14` rad. ON maximum
yaw, speed, torque, and source-signed power gaps are `0.001481` rad,
`0.000193` rad/s, `23.15` N m, and `0.1273` W, all within the existing
published-case limits. The ON source holds its heading coefficients for
small yaw changes while Python interpolates continuously; this numerical
difference remains explicit in the parity ledger.

Ordinary Python CI uses the two NPZ files and bundled HDF5. The existing
fresh MATLAB workflow continues to read newly generated source CSVs and
HDF5. These strengthen two already counted OSWEC cases.

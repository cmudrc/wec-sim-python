# OSWEC hinged body damping in an irregular sea

`source.npz` contains two local MATLAB R2025b runs of a derived, 150 s
`Nonhydro_Body` input. The first has zero flap body damping; the second uses
`diag([100000, 0, 0, 0, 2000000, 0])` in SI units. Both retain the zero PTO
damping, original hinge and body properties, and pinned `OSWEC.slx`. They use a
seed-1 PM sea with `Hs=2.5` m and `Tp=8` s, `simu.dt=0.025` s, 40 s radiation
memory, and the original 100 s ramp. Only these settings and the 150 s end
time differ from the published regular-wave input.

The MATLAB core revision is `0753b2e47f2457c078751dcfe5d251d1767b80ab`;
the Applications revision is `d53d4d4c9eda2581f04204f5d394a6ef84bb099e`.
Pinned `oswec.out` has SHA-256
`26339b72b21c9212b842e38db48c4150f60fbc43588a90e96b3c7daeb7937d99`.
Running its Applications `bemio.m` with the pinned core produces the source
HDF5 (SHA-256 `290fe562d5efe20eafb7b3e68a928eb0f7f40df190c721de6ed0e55eb1778564`).
The portable test uses the already bundled `oswec_multi_wave/oswec.h5`:
all 68 dataset paths match; four IRF arrays differ by at most `1.46e-11`
in their stored units, and the other datasets match exactly.

To regenerate, copy the pinned `Nonhydro_Body` case beside an
`_Common_Input_Files/OSWEC` directory containing that BEMIO HDF5 and the
pinned STL files. Make two copies of the case. In both inputs, set `endTime`
to 150, `dt` to 0.025, `cicEndTime` to 40, and change `waveClass('regular')`
to `waveClass('irregular')`; retain height 2.5 and period 8 and add direction
0, spread 1, spectrum type `PM`, and phase seed 1. Set flap
`linearDamping` to zero or the matrix above. Run
`tests/matlab/oswec_hinge_body_damping_pm_reference.m` on each copy; it checks
the settings and exports body, PTO, wave, phase, and damping CSVs. Stack the
control and damped CSV arrays in that order under their filenames without
the `.csv` suffix, using `numpy.savez_compressed`.

The packed fixture SHA-256 is
`cdba2e506c2b002863a4d8f07fd5d17c1f15b62dbae03d125cd7042722aa1722`.
Each arm has 6,001 time samples. The two source seas are identical, and the
damper changes source pitch by `0.052503` rad. Independent public `WEC.run`
matches the damper-induced surge, heave, and pitch changes within
`0.438` mm, `0.184` mm, and `0.0000944` rad; their speed-change errors are
`0.447` mm/s, `0.250` mm/s, and `0.0000948` rad/s. All pass the unchanged
differential limits from the regular-wave body-damping comparison. This
tests the radiation-memory path of a derived setting; `Nonhydro_Body` was
already counted, so the published-case badge does not change.

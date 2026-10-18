# RM3 End Stops internal solver times

The adjacent NPZ records a local R2025b run of the pinned WEC-Sim source
`0753b2e4` and Applications input `d53d4d4`, using the generated RM3 HDF5.
The published 400 s, 0.1 s `ode45` case was copied to a disposable directory.
Only Simulink time logging was enabled. The saved PTO trajectory is exactly
equal to `rm3_end_stops_maxstep_400s.npz`'s published arm at all 4,001 samples.

`tests/matlab/rm3_end_stops_maxstep_audit.m` reproduces the export. Its full
duration published arm saves the time log, both bodies' output accelerations,
applied added-mass forces, and post-shift added-mass matrices. The portable
test checks that every output time occurs in the solver log, then compares
the saved-grid delayed-acceleration reconstruction with the applied force.

The log has 81 off-grid times. Exactly 78 output samples have an added-mass
recurrence residual above 0.001 N, and each is within two output intervals
of an off-grid time. All other residuals are below `5e-6` N. This identifies
why a reconstruction using only saved 0.1 s accelerations fails after contact.
It does not establish an independent Python motion gate or change the
`End_Stops` qualified status.

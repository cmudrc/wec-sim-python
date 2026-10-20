# Published RM3 End Stops body trajectories

`rm3_end_stops_published_full.npz` contains the published 400 s, 0.1 s
`ode45` run from WEC-Sim `0753b2e47f2457c078751dcfe5d251d1767b80ab` and
Applications `d53d4d4c9eda2581f04204f5d394a6ef84bb099e` on local MATLAB
R2025b Update 7. The source input and Simulink model were copied to a
disposable directory; the only model change was enabling time logging.
`tests/matlab/rm3_end_stops_maxstep_audit.m` generated the published arm.

The fixture stores both bodies' surge, heave, and pitch positions and speeds,
the PTO stroke, speed, and internal force, wave elevation, and time. The new
PTO CSV is identical to the previously pinned published PTO array at all
4,001 samples. The 4,082 logged accepted solver times are identical to
`rm3_end_stops_internal_steps.npz`. Source CSV SHA-256 values:

- `published_body1.csv`: `d486b7d16c29902c16458f450162ffb54732c1cb82bdb49ccf624ca0bed23778`
- `published_body2.csv`: `0e79e4f3643eeb70c2d6d0e84182c04fbc3ba39ef9ce937e027b03cafd859f80`
- `published_pto.csv`: `683849829f7de3a6a92cd9c6a284cf66d8a817f94980a71033db0a145fd663b2`

MATLAB regenerated `rm3.h5` from the pinned WAMIT output. Its 84 dataset
paths match `examples/data/rm3_end_stops.h5`; active regular-wave
hydrodynamic datasets agree to saved precision. The files differ in
state-space radiation fits, which this regular-wave run does not use.

`test_rm3_end_stops_source_timing.py` compares independently advanced Python
motion against this fixture. The source-state ode45 controller independently
pairs both bodies' active positions and speeds, PTO motion and force, and work
through 120 s; it does not ingest source solver times. Its full-duration arm
supplies only MATLAB's accepted step times, not body states or forces. That
timing condition keeps the published case qualified until Python independently
selects the steps and passes the same full-duration gates.

A separate 90 s local source run reproduced the first 901 published PTO
samples and the first 903 accepted step times. MATLAB's Solver Profiler
reported zero zero-crossing events, zero solver exceptions, and one solver
reset at startup. Its first extra accepted steps occur at 84.380081548 and
84.444324664 s. The profiler export is a disposable diagnostic and is not
used as a test input.

# RM3 End_Stops MaxStep-only MATLAB source

`rm3_end_stops_maxstep_120s.npz` holds all 1,201 samples from the
[`maxstep_0025` arm](https://github.com/cmudrc/wec-sim-python/actions/runs/38849501711)
of the fresh R2025b audit at `cmudrc/wec-sim-python` commit `241e73add`.
The source uses WEC-Sim `0753b2e` and Applications `d53d4d4` with the
published 0.1 s input/output grid, wave, PTO, and 100 s ramp. Only the
source `ode45` MaxStep changes to 0.025 s; the run ends at 120 s.

The archived artifact CSV SHA-256 values are:

- `body1`: `20c577a5cededee805b30ca0ebeb4b031764a9bc2d902eacad3f91cefef4d1d1`
- `body2`: `b9199c593bb276409711d3a7ccc75bd4517d4f7d0bed2d25f9c17bafa3517ab4`
- `pto`: `0626a88d430c7574fc6321583c36c68bc42eb115c8c576beb6a9e7de46f315e5`
- `wave`: `2cd9061fbe8297007420c9c8460f80c5116b6353c2eb5f615aa5bc7c50e074f0`

The portable test advances Python from the HDF5 and case configuration; no
source position, velocity, or force history enters that solve. This is a
derived solver-setting comparison, so `End_Stops` remains qualified.

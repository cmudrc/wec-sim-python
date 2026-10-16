# Published RM3 End_Stops prefix

`rm3_end_stops_published_prefix.npz` contains the first 120.0 s, including
1,201 output samples, from the pinned R2025b MATLAB
[`RM3_END_STOPS` run](https://github.com/cmudrc/wec-sim-python/actions/runs/38367347607).
It stores time, each body's surge/heave/pitch position and velocity, and PTO
stroke/speed/force. The pinned source is WEC-Sim `0753b2e` with Applications
`d53d4d4`; the audit changes only the copied input's end time to 120 s.

The source CSV SHA-256 values before trimming were:

- `body1`: `b9d23dc098f7dfaec03014611067ce06c9bd1d0e7e45d6d1b954a792e1d01256`
- `body2`: `4c6fa3a6c90fa66bbee5efe65843641328d88a4cb2f9f78b257a58d8833072be`
- `pto1`: `9efad6333ca7c2d2e6a47e9328131067fe3bea360cd4e242a07f625e89afe11c`

The fixture is an identity check for the unchanged MATLAB arm. The separate
`rm3_end_stops_full.npz` fixture supplies the previously gated refined
source trajectory; neither is fed into an independent Python solve.

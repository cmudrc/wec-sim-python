# RM3 End_Stops full-duration MaxStep-only source

`rm3_end_stops_maxstep_400s.npz` contains the 4,001 saved samples from the
[fresh R2025b 400 s audit](https://github.com/cmudrc/wec-sim-python/actions/runs/38854937372)
at commit `b783e0d1d`. WEC-Sim `0753b2e47f2457c078751dcfe5d251d1767b80ab`
and Applications `d53d4d4c9eda2581f04204f5d394a6ef84bb099e` are pinned.
Both source arms keep the published 0.1 s input/output grid, regular wave,
PTO, and 100 s ramp. The changed arm sets only `ode45` MaxStep to 0.025 s.
Its first 120 s reproduces the separately archived changed-arm fixture;
the unchanged arm reproduces the published prefix. The two wave histories
match at every saved sample.

The source artifact CSV SHA-256 values are:

- Changed float: `e69eda4f1740b906dd6a1fb9f20318e4e628cc818b47d67e85ccf7901e477`
- Changed spar: `82bb56d39c26263a07d3c355e737631be19e188b269b11f738aa24fa76e17f5f`
- Changed PTO: `a58b00cae505bcc1730d5370344e5265d278533d7b01fc30eaef0394fd14b984`
- Published PTO: `683849829f7de3a6a92cd9c6a284cf66d8a817f94980a71033db0a145fd663b2`
- Wave (both arms): `a265719696773b356e697d02a5a79b2d88da855bd1a9086b12dbcb244ba9f119`

The compressed fixture SHA-256 is
`4a01cc6119c5ba1606e3ccdef3001cfc5b941a0c5881fc9a7631e8f7780f358d`.
The portable test advances the physical Python default from the HDF5 and
case configuration; MATLAB states and loads enter only the comparison. This
is a derived solver-setting gate. The unchanged published case remains
qualified.

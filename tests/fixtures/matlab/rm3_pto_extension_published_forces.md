# RM3 PTO Extension published-step forces

`rm3_pto_extension_published_forces.npz` records both bodies from each
unchanged 30 s, 0.1 s `ode4` input in the pinned MATLAB
`RM3_PTO_Extension/float` and `spar` applications. Columns are time, heave,
heave speed, heave acceleration, radiation-damping force, and added-mass
force. The HDF5 was regenerated through pinned BEMIO. The source revisions
are WEC-Sim `0753b2e47f2457c078751dcfe5d251d1767b80ab` and
Applications `d53d4d4c9eda2581f04204f5d394a6ef84bb099e`, with MATLAB
R2025b Update 7 on arm64 macOS. The export script is
[`tests/matlab/rm3_pto_extension_force_audit.m`](../../matlab/rm3_pto_extension_force_audit.m).

The fresh body trajectories match the archived published-step motion fixture
within `1e-13` m and `1e-14` m/s. The portable test reconstructs the active
heave radiation convolution from source velocity and the checked-in matching
RM3 IRF within `0.005` N; it reconstructs the delayed added-mass force from
the source acceleration history within `0.025` N. The largest logged
forces reach 3.095 MN and 13.885 MN for float and 17.319 kN and 2.178 MN
for spar, respectively. These checks isolate the major-step force laws.

Source motion and acceleration are inputs to these force reconstructions.
They do not establish independent published-step Python motion parity; the
float and spar cases remain qualified.

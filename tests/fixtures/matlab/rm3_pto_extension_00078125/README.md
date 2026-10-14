# RM3 PTO extension next-step audit

Local MATLAB R2025b Update 7 used pinned WEC-Sim core
`0753b2e47f2457c078751dcfe5d251d1767b80ab` and Applications
`d53d4d4c9eda2581f04204f5d394a6ef84bb099e`. A disposable copy of each
published `RM3_PTO_Extension` float and spar input changed only `simu.dt` from
the archived 0.003125 s to 0.0015625 s and 0.00078125 s. End time remained
30 s, the wave was `noWaveCIC`, the PTO stiffness and damping stayed zero,
and the extension target stayed 5 m. Both cases used the same common-RM3
HDF5 (SHA-256
`b8dd385f5db2a53fa3ffb971b105ea4fe8a2685ffd00bab0ae6050f366f69768`).
MATLAB exited successfully after writing four complete trajectories.

`source.npz` contains the nine source columns at every output step: time,
float heave/speed, spar heave/speed, PTO stroke/speed/force/power. Load each
CSV with `numpy.loadtxt(..., delimiter=",")` and save the arrays under
`float_0p0015625`, `float_0p00078125`, `spar_0p0015625`, and
`spar_0p00078125` using `numpy.savez_compressed`. Its SHA-256 is
`fc4307283892e4d5c18b26293c7039c5d543bdd0a26dcf9fb7b54fcf94c433b8`.
The source CSV SHA-256 hashes in that key order are:

| Source | SHA-256 |
| --- | --- |
| float 0.0015625 s | `2a4d771b8a9549f45a8744d76eb20fbec296be3bc7da6e180b532ed539e1ef25` |
| float 0.00078125 s | `c982f424cab25f59b4893f81c0dcfc1a790f008c7749908c711cc47861ee25d4` |
| spar 0.0015625 s | `55f65f7c720dca681ffd40a4799c123a214d22ce64c5995100c2d306b29b551b` |
| spar 0.00078125 s | `b1cd1073bc05dd0d77d1bef38aff79af85e09780de43e627ff50a8783ae7d18e` |

`gate_spec_pre_grid_correction.json` preserves the numerical screens written
before either run. Before evaluating a complete case, the comparison grid
was corrected from 0.003125 s to 0.00625 s because the oldest source history
has samples only every 0.00625 s. `gate_spec.json` records that indexing
correction; no limits or required decreasing trends changed. `metrics.json`
contains the resulting full-duration source self-changes, extrapolation
spreads, independent Python differences, and pass/fail results.

The Python test independently advances each public physical `WEC` model from
the same HDF5 with no saved source body motion or force as input. This is a
refined-step convergence diagnostic. The published 0.1 s inputs remain
qualified; these data do not promote their case badges.

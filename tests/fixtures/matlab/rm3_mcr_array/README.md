# RM3 MCR Option 1 array histories

`source.npz` is a lossless column selection from the successful [pinned
R2025b MATLAB run](https://github.com/cmudrc/wec-sim-python/actions/runs/37787926107)
of the published `RM3_MCROPT1` eight-condition array. That run used WEC-Sim
core `0753b2e47f2457c078751dcfe5d251d1767b80ab` and Applications
`d53d4d4c9eda2581f04204f5d394a6ef84bb099e`. The complete CSV exports
remain in its `matlab-RM3_MCR_ARRAY-baseline` artifact. The NPZ SHA-256 is
`2de50eb304b39c7fe37f08afc9d55168e9d461655ce3b29ae5bceb77f8a36ed1`.

An independent local MATLAB R2025b Update 7 run on 2026-10-15 used the same
pinned core and Applications revisions and regenerated `rm3.h5` with pinned
BEMIO. `reference_model_baseline('RM3_MCR_ARRAY')` completed all eight 400 s
conditions. Against this archive, the maximum absolute difference across the
selected body motion channels was `8.35e-13`; PTO stroke and speed differed by
at most `1.60e-14` m and `6.00e-15` m/s, PTO force and power by at most
`1.21e-8` N and W, and the eight-case mean-power summary by at most
`9.90e-10` W. These differences are at CSV serialization precision.

`summary` is the original 8×6 source table: wave height, period, PTO damping,
PTO stiffness, mean signed mechanical power, and repeated damping. For each
`case1` through `case8`, the 4,001×17 array selects these source CSV columns
without rounding:

| NPZ columns | Source |
| --- | --- |
| 0 | body 1 time |
| 1–6 | body 1 surge/heave/pitch position, then velocity |
| 7–12 | body 2 surge/heave/pitch position, then velocity |
| 13–16 | PTO stroke, speed, axial force, signed mechanical power |

The paired public `WEC.floating_joint` test regenerates each incident sea and
advances both bodies from the published parameters and bundled RM3 HDF5. It
passes no saved body motion or force to Python. The optional
`added_mass_scheme="simulink_delay"` reproduces the source feedback timing;
the ordinary implicit scheme remains the physical default. The fixture
deepens output gates for eight already counted MCR inputs and does not change
the case badge.

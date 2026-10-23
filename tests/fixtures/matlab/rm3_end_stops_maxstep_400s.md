# RM3 End_Stops corrected validation baseline

This is the **repository-owned corrected validation case** for RM3
`End_Stops`. Starting from the pinned Applications input, the
[in-repository MATLAB runner](../../matlab/rm3_end_stops_maxstep_audit.m)
copies `End_Stops` into a disposable directory and changes only the copied
Simulink model's `ode45` `MaxStep` from 0.1 to 0.025 s. It retains
`simu.dt=0.1` s, the 400 s duration, 100 s ramp, regular wave, PTO, and
0.1 s output grid. No WEC-Sim or Applications checkout is edited.

The 400 s `maxstep_0025` arm is the MATLAB baseline for
[`test_rm3_end_stops_maxstep_public.py`](../../test_rm3_end_stops_maxstep_public.py).
That test independently advances the public Python physical solver from the
pinned HDF5 and case settings and gates both active bodies' motion, PTO
stroke, speed, force, and ordinary damper work. It does not feed MATLAB
states, forces, events, or accepted step times into Python. The unchanged
0.1 s `published` arm is retained to measure source solver sensitivity. It
is not the corrected validation baseline and remains qualified in the
published-case inventory and badges.

To regenerate the baseline, use the pinned WEC-Sim and Applications revisions
below and run `rm3_end_stops_maxstep_audit(400)` as in the
[`rm3-end-stops-maxstep-full-audit.yml`](../../../.github/workflows/rm3-end-stops-maxstep-full-audit.yml)
workflow. Compare the generated CSVs with this record using
[`test_rm3_end_stops_maxstep_audit.py`](../../test_rm3_end_stops_maxstep_audit.py)
and then run the independent portable Python gate. The regeneration check
compares all 4,001 saved samples of both bodies, the PTO, and wave against
the checked fixture, allowing at most 1 mm for body translation and PTO
stroke, 1 mm/s for speed, 10 kN for PTO force, and `1e-12` m for wave.
The local R2025b rerun and its cross-platform checks are recorded below.

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
case configuration; MATLAB states and loads enter only the comparison.

## Local R2025b platform audit

On 2026-10-17, the same `rm3_end_stops_maxstep_audit.m` ran at both 120 and
400 s on macOS R2025b. The pinned WEC-Sim source was
`0753b2e47f2457c078751dcfe5d251d1767b80ab`; Applications was
`d53d4d4c9eda2581f04204f5d394a6ef84bb099e`. The End Stops model and
RM3 common inputs were copied into a disposable case tree before the run, so
the pinned checkouts were not changed. The predeclared cross-platform limits
were 1 mm for active position and PTO stroke, 1 mm/s for speed, 10 kN for PTO
force, and `1e-12` m for the wave. All 120 and 400 s comparisons passed.

| 400 s comparison with archived Linux R2025b | Maximum absolute difference |
| --- | ---: |
| Published PTO stroke, speed, and force | 0 (byte-identical CSV) |
| MaxStep-only active body position | `1.03e-13` m or rad |
| MaxStep-only active body speed | `9.99e-15` m/s or rad/s |
| MaxStep-only PTO stroke | `2.05e-15` m |
| MaxStep-only PTO speed | `1.30e-14` m/s |
| MaxStep-only PTO force | `1.60e-7` N |
| Wave elevation | 0 (byte-identical CSV) |

The local 400 s published-versus-MaxStep-only gap reaches `0.2854494` m PTO
stroke, `0.8601823` m/s PTO speed, and `8.006595` MN PTO force. The source
platform therefore does not account for the previously measured published
step sensitivity. No Python motion gate or badge count changes from this
source-to-source check.

Local CSV SHA-256 values: published PTO
`683849829f7de3a6a92cd9c6a284cf66d8a817f94980a71033db0a145fd663b2`,
MaxStep-only float
`fa582914862208fb3ea7fa2f146ab9bfaa9c0e27f8fa609424f1f8d823a07ea9`,
spar `552fbf305ef8e74e24c8dc74a2688abf2b97e8b117b33cf41f22c92e8241c560`,
PTO `61234a9237b6ab0fb707943338369311cdb3e6fb5a3037bf868c3a2b18770716`,
and wave `a265719696773b356e697d02a5a79b2d88da855bd1a9086b12dbcb244ba9f119`.

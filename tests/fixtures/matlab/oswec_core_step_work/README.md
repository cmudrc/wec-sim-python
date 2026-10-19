# Core OSWEC fixed-step PTO-work audit

The published `examples/OSWEC` input uses a 0.1 s fixed `ode4` step for
400 s. A new 1% accumulated PTO-work comparison against that published trace
failed: independent public Python produced 18,677.37 J versus MATLAB's
18,885.16 J, a 1.1003% gap. The existing full-duration motion, PTO force,
and instantaneous-power gates still pass. This audit changes `simu.dt`,
which also changes WEC-Sim's force and output clocks. A separate
[solver-only audit](../oswec_core_fixedstep_only/README.md) holds those clocks
at 0.1 s. Neither audit changes the strict case count.

The local source is WEC-Sim core `0753b2e47f2457c078751dcfe5d251d1767b80ab`
under R2025b Update 7 on Mali. Its pinned OSWEC HDF5 and the bundled Python
copy both have SHA-256 `28aa028d404a88f877924bc435e7ad3a4c50e26f9706c1f39ef7324e5ec60750`.
The unedited `OSWEC.slx` and `wecSimInputFile.m` hashes are
`1d5f16c2a56337d654b727cc5f828996ebdfb8cb3d23e671cb0013334e52c26a`
and `13202805ad8b259067b90f3b9fffd1d0918d69f38abb08e07da25b81a123bf08`.
The [source script](../../../matlab/oswec_core_step_work_audit.m) copies the
case to a scratch directory, fixes `waves.phaseSeed=1` to reproduce the
published fixture, disables the post-run plots and Explorer, and changes
only `simu.dt` to 0.05, 0.025, or 0.0125 s. Each arm ran in its own MATLAB
process. The 0.025 s control exactly reproduced the archived PTO, wave, and
flap arrays. Two separate 0.0125 s runs produced byte-identical CSV files;
the checked-in source script reproduces that arm separately.

| Full-run absorbed work | MATLAB | Change from prior MATLAB step |
| --- | ---: | ---: |
| Published 0.1 s | 18,885.16 J | — |
| 0.05 s | 18,758.89 J | −0.6686% |
| 0.025 s | 18,730.10 J | −0.1535% |
| 0.0125 s | 18,723.17 J | −0.0370% |
| Independent Python at 0.025 s | 18,718.19 J | −0.0636% from 0.025 s MATLAB |
| Independent Python at 0.0125 s | 18,720.22 J | −0.0130% from 0.0125 s MATLAB |

The last three MATLAB steps have a successive work-difference ratio of
4.153 (observed order about 2.05). Independent public Python at 0.025,
0.0125, and 0.00625 s has a ratio of 4.0003. Second-order extrapolation gives
18,720.856 J for MATLAB and 18,720.903 J for Python, a 0.000252% gap.
Across the full run and three ramp/steady intervals, the new portable
[convergence regression](../../../test_oswec_core_eighth_convergence.py)
checks the incident sea and passive PTO law, full-motion errors at 0.0125 s,
the unchanged 0.1% refined-step work limit, ratios between 3.5 and 4.5,
and a 0.01% extrapolated-work limit. These are empirical finite-step results,
not an exact convergence proof or a passed published 0.1 s work gate.

All 500 phase-bearing wave components come from the published fixture.
MATLAB wave elevation is exactly identical at all shared 0.1 s samples
across the three source steps. The PTO force and power of each refined trace
obey the unchanged 12,000 N m s/rad passive-damper law within `1e-8` in
their units. The 0.1-to-0.05 s source change moves PTO speed by up to
0.002742 rad/s and torque by 32.90 N m; the 0.05-to-0.025 s changes fall to
0.000636 rad/s and 7.64 N m. The corresponding peak flap-pitch changes are
0.001573 and 0.000401 rad.

The [portable gate](../../../test_oswec_core_step_work_audit.py) advances
Python independently for all 16,001 quarter-step samples from the pinned
HDF5 and wave phases. Its peak surge, heave, and pitch position differences
from 0.025 s MATLAB are 0.827 mm, 0.177 mm, and 0.168 mrad. PTO speed,
torque, and instantaneous absorbed-power differences peak at
0.000232 rad/s, 2.785 N m, and 0.969 W. The full-run and three ramp/steady
interval work gaps are 0.064%, 0.064%, 0.062%, and 0.066%, all below the
0.1% CI gate. This is a derived-step physical-motion gate; it does not
retroactively pass the failed published-step 1% work limit or prove MATLAB
step convergence.

`eighth_source.npz` (SHA-256
`33b6455b79162c1f161f63cab23a290556a9ed8fc0d77ddab5c89140090671b2`)
holds the 0.0125 s PTO, incident wave, and active flap histories.
`source.npz` (SHA-256
`4558bc3b9485a4a19cdabf109239af96c8111a351eea0e14535e0ff85013490f`)
holds the 0.05 and 0.025 s source PTO, wave, and active flap histories.
The original 0.1 s trace is the existing `oswec_irregular_wave.npz` fixture.

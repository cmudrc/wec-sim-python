# Core OSWEC fixed-step PTO-work audit

The published `examples/OSWEC` input uses a 0.1 s fixed `ode4` step for
400 s. A new 1% accumulated PTO-work comparison against that published trace
failed: independent public Python produced 18,677.37 J versus MATLAB's
18,885.16 J, a 1.1003% gap. The existing full-duration motion, PTO force,
and instantaneous-power gates still pass. This audit checks whether changing
only MATLAB's solver step accounts for the work difference; it does not change
the strict case count.

The local source is WEC-Sim core `0753b2e47f2457c078751dcfe5d251d1767b80ab`
under R2025b Update 7 on Mali. Its pinned OSWEC HDF5 and the bundled Python
copy both have SHA-256 `28aa028d404a88f877924bc435e7ad3a4c50e26f9706c1f39ef7324e5ec60750`.
The unedited `OSWEC.slx` and `wecSimInputFile.m` hashes are
`1d5f16c2a56337d654b727cc5f828996ebdfb8cb3d23e671cb0013334e52c26a`
and `13202805ad8b259067b90f3b9fffd1d0918d69f38abb08e07da25b81a123bf08`.
The [source script](../../../matlab/oswec_core_step_work_audit.m) copies the
case to a scratch directory, fixes `waves.phaseSeed=1` to reproduce the
published fixture, disables the post-run plots and Explorer, and changes
only `simu.dt` to 0.05 or 0.025 s. Each arm ran in its own MATLAB process.
Both refined outputs were independently reproduced byte for byte using the
source script after the original local runs.

| Full-run absorbed work | MATLAB | Change from prior MATLAB step |
| --- | ---: | ---: |
| Published 0.1 s | 18,885.16 J | — |
| 0.05 s | 18,758.89 J | −0.6686% |
| 0.025 s | 18,730.10 J | −0.1535% |
| Independent Python at 0.025 s | 18,718.19 J | −0.0636% from 0.025 s MATLAB |

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

`source.npz` (SHA-256
`4558bc3b9485a4a19cdabf109239af96c8111a351eea0e14535e0ff85013490f`)
holds the 0.05 and 0.025 s source PTO, wave, and active flap histories.
The original 0.1 s trace is the existing `oswec_irregular_wave.npz` fixture.

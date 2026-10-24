# RM3 PTO Extension corrected validation baseline

This repository adopts separate corrected validation baselines for the
published `RM3_PTO_Extension/float` and `RM3_PTO_Extension/spar` 30 s no-wave
free decays. Starting from WEC-Sim core
`0753b2e47f2457c078751dcfe5d251d1767b80ab` and Applications
`d53d4d4c9eda2581f04204f5d394a6ef84bb099e`, the
[MATLAB runner](../../../matlab/rm3_pto_extension_corrected.m) copies each
published model and input into a disposable folder. It changes only the
copied input's `simu.dt` from 0.1 to 0.00078125 s. Duration, no-wave
convolution mode, hydrodynamic data, zero PTO stiffness/damping, and 5 m
extension target stay fixed. No upstream checkout or repository is edited.

The original 0.1 s cases remain qualified because their physical-default
Python trajectories miss the published motion gates. This copied-input
baseline follows the repository's adopted corrected-validation policy in
`PARITY.md`; it is not a claim that the pinned upstream application was
changed. The historical published-step source-ode4 and force-law gates remain
intact.

`source.npz` contains the original 301 reporting samples, without
interpolation, from each full 38,401-sample 0.00078125 s MATLAB run. Its
`float` and `spar` arrays have columns time, float heave/speed, spar
heave/speed, PTO stroke/speed/force/power. The dense source histories remain
checked in the prior
[`rm3_pto_extension_00078125/source.npz`](../rm3_pto_extension_00078125/source.npz)
fixture. Both runs were generated locally with MATLAB R2025b Update 7 and
the pinned common-RM3 HDF5 recorded there. The compact corrected archive
SHA-256 is
`13bee7428cd48c79b507c93a10a77caa3c3f456606d735c193504886bd0bb235`;
the dense source archive SHA-256 is
`fc4307283892e4d5c18b26293c7039c5d543bdd0a26dcf9fb7b54fcf94c433b8`.

The [portable gate](../../../test_rm3_pto_extension_corrected.py) first checks
the complete source self-change between 0.0015625 and 0.00078125 s. The
predeclared limits are 1 mm heave and 3 mm/s speed for float, 0.3 mm and
0.2 mm/s for spar. The measured changes are 0.610 mm / 2.449 mm/s and
0.125 mm / 0.0764 mm/s. It then independently advances each public physical
Python `WEC` for the same 30 s at 0.00078125 s from the pinned HDF5 and
published 5 m initial extension, without saved MATLAB states, accepted
times, force, or PTO history as inputs. Direct active-body heave/speed and
PTO stroke/speed must satisfy those same predeclared limits at every fine
sample; inactive body heave/speed must stay within 1e-8 m and m/s, and zero
PTO force/power within 1e-7 N and W. The gate also compares the Python
trajectory with this distinct 301-sample corrected archive and preserves the
earlier source-extrapolation checks. Python's other body axes keep the
published-case stationary-axis limits. The hosted regeneration additionally
compares all six position and speed axes of each independently advanced body
against the fresh complete MATLAB body histories.

The [hosted MATLAB workflow](../../../../.github/workflows/rm3-pto-extension-corrected.yml)
regenerates both source steps for each copied input on R2025b, requires every
new source sample to reproduce the checked archive within 1e-6, verifies the
generated HDF5 active mass/restoring/radiation data against the bundled
physical input, and runs the same independent Python gate against the fresh
source. This regeneration is required before promoting either case badge.

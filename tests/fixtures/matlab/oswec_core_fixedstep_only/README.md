# Core OSWEC solver-only fixed-step audit

The pinned `examples/OSWEC` application runs `ode4` at 0.1 s for 400 s.
The [source script](../../../matlab/oswec_core_fixedstep_only.m) copies the
case, fixes the same phase seed used by the published fixture, and runs three
local R2025b Update 7 arms. It changes **only** Simulink `FixedStep` after
WEC-Sim initialization: 0.1, 0.05, or 0.025 s. `simu.dt`, `dtOut`, `cicDt`, and
`nonlinearDt` remain 0.1 s in every arm. The source checkout is core commit
`0753b2e47f2457c078751dcfe5d251d1767b80ab` on Mali.

Invoke the script in MATLAB's base workspace, with a fresh output directory
for each arm:

```matlab
addpath('/path/to/wec-sim-python/tests/matlab');
coreRoot = '/path/to/wec-sim-matlab-reference';
outputDir = '/tmp/oswec-fixedstep-0025';
fixedStep = 0.025;
oswec_core_fixedstep_only
```

The fresh 0.1 s run reproduces the archived published fixture: wave samples
are identical, maximum flap-pitch difference is `1.03e-14` rad, and integrated
PTO work differs by `1.5e-10` J. All three runs have identical archived wave
elevation at every 0.1 s output sample and obey the unchanged 12,000 N m
s/rad linear PTO law. A separate 0.05 s rerun using the checked-in script
reproduced its PTO, wave, body, and settings CSV arrays exactly.

| MATLAB solver step; WEC-Sim clocks all 0.1 s | Full-run absorbed work | Change from prior step |
| --- | ---: | ---: |
| 0.1 s, published | 18,885.16 J | — |
| 0.05 s | 18,700.86 J | −0.9759% |
| 0.025 s | 18,661.34 J | −0.2113% |

With its own 0.1 s physical integration, Python gives 18,677.37 J: 1.1003%
below the published MATLAB work, but 0.0859% above the 0.025 s solver-only
source. Its four ramp/steady cumulative-work gaps to the latter are at most
0.1612%; full-duration maximum surge, heave, and pitch gaps are 21.3 mm,
3.80 mm, and 4.27 mrad. The separate all-clock refinement gives 18,730.10 J
at `simu.dt=0.025` s, showing that changing WEC-Sim's clocks partly offsets
the solver-only effect.

The [portable gate](../../../test_oswec_core_fixedstep_only.py) checks the
source clocks, waves, PTO law, fresh baseline identity, source self-change,
and an independently advanced public Python trajectory and work against the
solver-refined source. This establishes the leading numerical mechanism for
the PTO-work difference. The 1% **published-step** work gate still fails;
the source is step-sensitive, and neither convergence nor a new strict case
gate is claimed. `source.npz` has SHA-256
`c6f047eddc3c566951e744bd40a67be6304fe341ac680b33d9850f7db29613ca`.

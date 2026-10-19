# OSWEC Wave Markers PTO-work step audit

The published 400 s `Wave_Markers/OSWEC` run is archived in
[`oswec_wave_markers.npz`](../oswec_wave_markers.npz) from [R2025b run
37932186068](https://github.com/cmudrc/wec-sim-python/actions/runs/37932186068).
That run pins WEC-Sim core `0753b2e47f2457c078751dcfe5d251d1767b80ab`
and Applications `d53d4d4c9eda2581f04204f5d394a6ef84bb099e`.

The local R2025b audit on Mali copied that application and its exact generated
`oswec.h5` (SHA-256
`79d9ab232958f3b3f9ec1bcdcd9a56c9fe8eb54f22e9de2714403860b8a51b32`)
to a disposable directory. It changed only `simu.dt` from 0.1 to 0.05 s and
set `simu.dtOut=0.1` s, preserving the published 400 s duration, 100 s ramp,
2.5 m/8 s regular sea, and 12,000 N m s/rad PTO damping. The [MATLAB audit
script](../../../matlab/oswec_wave_markers_step_audit.m) recreates the run.
Its source result has 4,001 shared PTO/body samples, 8,001 wave samples,
and wave elevation exactly equal to the published run at all shared times.
A second local execution of the checked-in audit script reproduced the PTO,
flap, and wave CSVs byte for byte (SHA-256 `9af3cf86cdfcf948b74c350fa715872a2b0e9556dc24d55dd047365c5f3884db`,
`5312b2fe0d65a0087069a306277cca183f89d5db9f5574444eb66d3b1e435753`,
and `9640f1d7c1fde43bfe08735c026a6dba77fce48141fa45ea98e484241629972b`).

`half_source.npz` retains the shared 4,001 time, wave, flap pitch/speed, and
PTO speed/force/source-signed-power samples without rounding. Its SHA-256 is
`31969cc668a6332e222530a3db81b47f6122c17586e935f6a03a659d5d03fe35`.
The portable [Python gate](../../../test_oswec_wave_markers_work_step.py)
rebuilds the full case from the bundled HDF5 and compares each interval.

| Interval | Published MATLAB work (J) | Half-step MATLAB work (J) | Python work (J) | MATLAB self-change | Python/half-step gap |
| --- | ---: | ---: | ---: | ---: | ---: |
| 0–100 s | 18,247.04 | 17,982.95 | 17,898.17 | −1.447% | −0.471% |
| 100–250 s | 75,040.68 | 73,815.96 | 73,426.00 | −1.632% | −0.528% |
| 250–400 s | 75,099.61 | 73,868.71 | 73,476.46 | −1.639% | −0.531% |
| 0–400 s | 168,387.33 | 165,667.63 | 164,800.62 | −1.615% | −0.523% |

The originally attempted 1% published-work limit fails by 2.130% overall.
The same 1% limit passes against the derived half-step source in every window.
One step refinement does not establish convergence, so the published energy
gap stays qualified and does not change the case badge.

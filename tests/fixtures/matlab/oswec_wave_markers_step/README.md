# OSWEC Wave Markers PTO-work step audit

The published 400 s `Wave_Markers/OSWEC` run is archived in
[`oswec_wave_markers.npz`](../oswec_wave_markers.npz) from [R2025b run
37932186068](https://github.com/cmudrc/wec-sim-python/actions/runs/37932186068).
That run pins WEC-Sim core `0753b2e47f2457c078751dcfe5d251d1767b80ab`
and Applications `d53d4d4c9eda2581f04204f5d394a6ef84bb099e`.

The local R2025b audit on Mali copied that application and its exact generated
`oswec.h5` (SHA-256
`79d9ab232958f3b3f9ec1bcdcd9a56c9fe8eb54f22e9de2714403860b8a51b32`)
to disposable directories. It changed only `simu.dt` from 0.1 to 0.05 and
0.025 s and
set `simu.dtOut=0.1` s, preserving the published 400 s duration, 100 s ramp,
2.5 m/8 s regular sea, and 12,000 N m s/rad PTO damping. The [MATLAB audit
script](../../../matlab/oswec_wave_markers_step_audit.m) recreates either run.
The results have 4,001 shared PTO/body samples and 8,001 or 16,001 wave samples,
and wave elevation exactly equal to the published run at all shared times.
A second local execution of the checked-in audit script reproduced the PTO,
flap, and wave CSVs byte for byte (SHA-256 `9af3cf86cdfcf948b74c350fa715872a2b0e9556dc24d55dd047365c5f3884db`,
`5312b2fe0d65a0087069a306277cca183f89d5db9f5574444eb66d3b1e435753`,
and `9640f1d7c1fde43bfe08735c026a6dba77fce48141fa45ea98e484241629972b`).
The same checked-in script reproduced the 0.025 s PTO, flap, and wave CSVs
byte for byte (SHA-256 `2cef84d1743a229439ec053a66d34ebd0e65c9bde6d902bc13a7a47c980aeedf`,
`015aa4b3919854c78038a183eadb7f52d85f018e781cabca5fa8bfd24aa7dfed`,
and `ffd59a675dc2e542e05f3ad41db23954f752ee045b6ec71c9aeadb632ef17b98`).

`half_source.npz` and `quarter_source.npz` retain the shared 4,001 time, wave,
flap pitch/speed, and PTO speed/force/source-signed-power samples without
rounding. Their SHA-256 hashes are
`31969cc668a6332e222530a3db81b47f6122c17586e935f6a03a659d5d03fe35`
and `86f2e811ab53c1453ee2baaebec71a3312acd50586cf4fb66f10996785c1755b`.
The portable [Python gate](../../../test_oswec_wave_markers_work_step.py)
rebuilds the full case from the bundled HDF5 and compares each interval.

| Interval | Published MATLAB (J) | 0.05 s MATLAB (J) | 0.025 s MATLAB (J) | Python (J) | MATLAB 0.1→0.05 | MATLAB 0.05→0.025 | Python/0.025 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0–100 s | 18,247.04 | 17,982.95 | 17,919.19 | 17,898.17 | −1.447% | −0.355% | −0.117% |
| 100–250 s | 75,040.68 | 73,815.96 | 73,522.47 | 73,426.00 | −1.632% | −0.398% | −0.131% |
| 250–400 s | 75,099.61 | 73,868.71 | 73,573.53 | 73,476.46 | −1.639% | −0.400% | −0.132% |
| 0–400 s | 168,387.33 | 165,667.63 | 165,015.19 | 164,800.62 | −1.615% | −0.394% | −0.130% |

The originally attempted 1% published-work limit fails by 2.130% overall.
The same 1% limit passes against both derived refined sources in every window.
Two step refinements show a decreasing MATLAB self-change and a shrinking
Python/source gap; they do not establish convergence or pass the published
energy limit. The gap stays qualified and does not change the case badge.

# OSWEC Wave Markers PTO-work step audit

The published 400 s `Wave_Markers/OSWEC` run is archived in
[`oswec_wave_markers.npz`](../oswec_wave_markers.npz) from [R2025b run
37932186068](https://github.com/cmudrc/wec-sim-python/actions/runs/37932186068).
That run pins WEC-Sim core `0753b2e47f2457c078751dcfe5d251d1767b80ab`
and Applications `d53d4d4c9eda2581f04204f5d394a6ef84bb099e`.

The first local R2025b audit on Mali copied that application and its exact generated
`oswec.h5` (SHA-256
`79d9ab232958f3b3f9ec1bcdcd9a56c9fe8eb54f22e9de2714403860b8a51b32`)
to disposable directories. It changed only `simu.dt` from 0.1 to 0.05 and
0.025 s and
set `simu.dtOut=0.1` s, preserving the published 400 s duration, 100 s ramp,
2.5 m/8 s regular sea, and 12,000 N m s/rad PTO damping. The [MATLAB audit
script](../../../matlab/oswec_wave_markers_step_audit.m) recreates these runs.
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

For the new 0.0125 s arm, local R2025b ran the pinned application's
`_Common_Input_Files/OSWEC/hydroData/bemio.m` on its pinned `oswec.out` in a
disposable directory. The regenerated HDF5 SHA-256 is
`096ea5e696cfdf8205ab8a24a2226c2e02669d3345d97fc5483d85cb3ec3edac`;
its file bytes differ from the earlier source HDF5. A full 400 s control
rerun at 0.05 s using this HDF5 reproduced every archived shared PTO, flap,
and wave sample **exactly**. The same script then ran 0.0125 s with all other
inputs and the 0.1 s output grid fixed. Its PTO, flap, and 32,001-row wave
CSV SHA-256 values are
`68be37200fc066459bd30b90f245d18e40a2170e3eb325d1896d194e9208e86c`,
`df021635cd0abdb91ef2451db320a1c8db445d2598f24f6cd44d253b1efec741`,
and `e4b1db63c6ed64f046d9e9d0650a043cbca20525da11b85cc6ac6aa49767c5c1`.
A second local execution of the checked-in audit script reproduced all three
CSV hashes exactly. The eighth-step wave is identical at all 4,001 published
output times.

`half_source.npz`, `quarter_source.npz`, and `eighth_source.npz` retain the shared 4,001 time, wave,
flap pitch/speed, and PTO speed/force/source-signed-power samples without
rounding. Their SHA-256 hashes are
`31969cc668a6332e222530a3db81b47f6122c17586e935f6a03a659d5d03fe35`
and `86f2e811ab53c1453ee2baaebec71a3312acd50586cf4fb66f10996785c1755b`,
and `db4cc9cdbcbe15f7aa9542373e4b00108a33e3c6f8522f46a577a1ad107eba18`.
The portable [Python gate](../../../test_oswec_wave_markers_work_step.py)
rebuilds the full case from the [exact source-generated HDF5](../oswec_wave_markers/README.md)
and compares each interval. The earlier audit used the bundled 10-heading
example HDF5; the corrected comparison uses the published 36-heading input.

| Interval | Published MATLAB (J) | 0.05 s MATLAB (J) | 0.025 s MATLAB (J) | 0.0125 s MATLAB (J) | Python (J) | Observed order | Python/second-order extrapolation |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0–100 s | 18,247.04 | 17,982.95 | 17,919.19 | 17,903.47 | 17,898.25 | 2.019 | +0.000176% |
| 100–250 s | 75,040.68 | 73,815.96 | 73,522.47 | 73,450.27 | 73,426.37 | 2.023 | +0.000237% |
| 250–400 s | 75,099.61 | 73,868.71 | 73,573.53 | 73,500.88 | 73,476.84 | 2.023 | +0.000231% |
| 0–400 s | 168,387.33 | 165,667.63 | 165,015.19 | 164,854.62 | 164,801.47 | 2.023 | +0.000228% |

The originally attempted 1% published-work limit still fails by 2.130% overall.
The same 1% limit passes against all three derived refined sources in every
window. The final source halving changes full-duration work by −0.0973%; its
observed second-order extrapolation is 164,801.09 J versus independent
Python's 164,801.47 J. The regression bounds on the observed step ratio
and extrapolated agreement codify these measured results; they were not
predeclared published-case parity limits. The finite-step trend does not prove
mathematical convergence or pass the published 0.1 s energy limit. The gap
stays qualified and the case badge does not change.

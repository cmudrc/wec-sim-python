# Published WaveBot adaptive control case

These inputs and outputs come from the successful [100 s R2025b source
run](https://github.com/cmudrc/wec-sim-python/actions/runs/37863811457),
`matlab-WAVEBOT_CONTROL-baseline`, at repo commit
`a487c8aeca1188ba7295237855e60197c0356896`. The workflow pins WEC-Sim
core `0753b2e47f2457c078751dcfe5d251d1767b80ab` and Applications
`d53d4d4c9eda2581f04204f5d394a6ef84bb099e`.

`waveBotBuoy.h5` is the source run's exact generated hydrodynamic file.
`Z_included.mat` and `f_vec.mat` are byte-for-byte copies of the two published
Applications `Load_Mitigating_Controls/ControlTests` inputs at the pinned
revision. Their SHA-256 hashes are, respectively,
`5962bbd5dd65abfc3e0293865ac86266ccd4a582a8da5c7b38c01b41b00074d1`,
`779f9ff94a05ff2720e14feecf5c38d7076dbb9e3c3c1bcfb3926ca50834f4f2`,
and `51ceb0118e8edd69d47a9812a039870aa5ab5e957c773ce85fcf78ea92d88c93`.

`published_100s.npz` retains float64 source samples without rounding. It
contains all 500 wave components; 10,001 samples of wave elevation, active
body position and speed, and six controller gains; the raw actuation time
series; and all 13 complex excitation estimates. Its SHA-256 is
`7d8d02db4989198e150152b1b6438751446b582ce63301f542d86be0c2528ebd`.
The source CSV hashes and extraction columns are:

| Source CSV suffix | SHA-256 | Archive data |
| --- | --- | --- |
| `components` | `ade4ed4553c7b4528a22f72ec14213dccec7060fe7ea18da79dcf698cb81bec4` | all columns |
| `wave` | `30c34f7c52f7267474a01440fa8b3c0813cbffc578fd218fcf63370c48110652` | all columns |
| `ControlTests_body1` | `b738257173d5d8d34f7a44c58eb65529818f58a4cd323baa00704e2e2faafc44` | time 0, position 1/3/5, speed 7/9/11 |
| `gains_real` | `1616037f926039d84ddc2ca2b27ff8bfe88a51eeef86b2d8bfa494a7c765fc89` | columns 1–6 |
| `F_actuation_real` | `1086778dd00274a1a325dc87ce4f1abec37572ebdf4476e5dfe8af0cc3914085` | all columns |
| `FexcOut_real` | `bb80322fcac68ee954662479c0fc7516c789fa41a4a27f52392b6fd771c8b1f8` | columns 1–615 |
| `FexcOut_imag` | `55cdea079862e695bfccc80e90b074555f13e4615fcae890e6a413b79e0860f2` | columns 1–615 |

Each row of the last two CSVs combines into a complex `(205, 3)` array in
Fortran order. `test_wavebot_power_control_portable.py` uses only the saved
wave phases as dynamic input. It independently generates the sea, advances
the 100 s plant and PI controller, estimates excitation, and optimizes all
13 gain updates. Its existing source-comparison limits are unchanged.
The published input's indexed linear damping is selected explicitly; the
physical dissipative default remains separate. This published case was
already counted in the WaveBot parity badge before this portable gate.

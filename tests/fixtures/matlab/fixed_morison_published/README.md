# Published fixed Morison monopile

`source.npz` packs the unrounded float64 columns of six CSVs from the
[pinned R2025b MATLAB run](https://github.com/cmudrc/wec-sim-python/actions/runs/37768083149)
(`matlab-MORISON_FIXED-baseline`, commit
`bde2332db1313f50a51318436dd40b0cd94e2eec`). It contains all 40,001
samples of the published 400 s, three-heading PM case. The archive SHA-256 is
`cc6e2579c0aeddaf4925940334cc135e4f093b708b59c7a48084d7cdcf171de4`.

| Archive key | Source CSV | Source SHA-256 |
| --- | --- | --- |
| `components` | `MORISON_FIXED_components.csv` | `d1da6b9f5df48e17eec4e625a55eba7469d86fbcec85ee1ad0fbded1dc227343` |
| `directions` | `MORISON_FIXED_directions.csv` | `fb581ce2fe25dcc31394e43c243ff2281627613fe0df5a543f712d4e9343cebf` |
| `wave` | `MORISON_FIXED_wave.csv` | `7223cb256ac91952441e90879e45bb7bfeedf3bcac42e65eb170708581dd48e9` |
| `force` | `MORISON_FIXED_morison_force.csv` | `e2915f949e9a0270c312614b9c420ef50486a37a21d5cdc4a2450d6fc831c057` |
| `monopile` | `MORISON_FIXED_morisonElement_body1.csv` | `a663c6ab7a62c99d39745d6b6d5c86a8180073cf8c00572b18bc2d7350e98c9a` |
| `tower` | `MORISON_FIXED_morisonElement_body2.csv` | `25aa640398864c7800b2a05e612aee9289161a1568e52bb4429e7987a5a152ad` |

`tests/test_morison_fixed_parity.py` regenerates the sea and force from the
published input and checks the full source records through `WEC.run`. Saved
source force and body motion are comparison targets, never simulation inputs.
Ordinary Python CI uses this archive when fresh MATLAB CSVs are unavailable;
the same test can use fresh CSVs through `WEC_SIM_MATLAB_MODEL_OUTPUT_DIR`.
The numerical limits are unchanged from the original paired test. This makes
an already counted case durable; it does not increase the case badge.

# RM3 direct linear generator source window

`rm3.h5` and `published_20s.npz` come from the successful [R2025b published
RM3_DD_PTO run](https://github.com/cmudrc/wec-sim-python/actions/runs/38507893734)
(`matlab-RM3_DD_PTO-baseline`, run commit
`d65d29c5eb88fbc02b891b2525223e24924f8656`). The workflow pins WEC-Sim
core `0753b2e47f2457c078751dcfe5d251d1767b80ab` and Applications
`d53d4d4c9eda2581f04204f5d394a6ef84bb099e`. The HDF5 is copied without
alteration (SHA-256 `0e2ba661c8851ac7885f17c36e47d30438dce9e5c8733f26c5baed4d6bd4cda5`).

The archive retains source samples 0–2000 (0–20 s, 0.01 s output) without
rounding or interpolation. Load the artifact's `body1`, `body2`, `pto1`, and
`drive` CSVs with `numpy.loadtxt(..., delimiter=",", max_rows=2001)` and
write these columns with `numpy.savez_compressed`:

| Archive key | Source columns |
| --- | --- |
| `time` | drive 0 |
| `body_heave`, `body_speed` | body1/body2 3 and 9, stacked as float/spar |
| `pto_stroke`, `pto_speed` | PTO 3 and 9 |
| `force`, `absorbed_power`, `friction_force` | drive 2, 1, and 3 |
| `phase_current`, `phase_voltage`, `electrical_power` | drive 4–6, 7–9, and 11 |

The source CSV SHA-256 hashes, in body1/body2/PTO/drive order, are
`5184a22120fa87e48530dd801cf3a74181887f8438c8f16427851340b135deee`,
`f631d17f1439548c5f0ff836996efb3e1fa52ec86bd7dfff56466f3f2e22b182`,
`a336c8af70ea2174491350280574553bba80195dd018f366257f1b797d24624a`,
and `2c1880820f10e37fda10b8164f6416cdbef4d6e3f2ac9947b7c204604d7bc5a5`.
The NPZ SHA-256 is
`412af0cfb1054583aec2d56147dc3dbcd600660ab228eeeb80eeab930402fea9`.

The portable test independently advances the public `WEC` body-local PTO and
three-phase generator at the published 0.0005 s integration step. It checks
both body heaves and speeds, PTO stroke and speed, force, absorbed mechanical
power, friction, electrical power, and three-phase current and voltage against
the source. The separate MATLAB workflow still gates the full 400 s case;
this 20 s fixture keeps a significant early dynamic window active in ordinary
Python CI without adding an 800,000-step run to every commit.

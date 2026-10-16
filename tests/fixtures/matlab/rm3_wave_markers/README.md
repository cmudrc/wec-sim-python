# RM3 published Wave Markers fixture

The pinned [R2025b published run](https://github.com/cmudrc/wec-sim-python/actions/runs/37934252549)
produced `matlab-WAVE_MARKERS_RM3-baseline` at run commit
`4f62dca899baa96505f557362350800b46bf1fa0`. The run uses the pinned
WEC-Sim core and Applications inputs recorded in `tests/reference_cases.csv`.
`source.npz` keeps all 4,001 times, the 25 world-XY locations and elevation
histories, origin elevation, and source wave parameters. It copies the CSV
values as float64 without rounding. No source motion or force enters the
portable marker test; the published RM3 body and PTO motion have a separate
portable gate in `tests/test_rm3_portable_parity.py`.

| Source CSV | SHA-256 |
| --- | --- |
| `WAVE_MARKERS_RM3_locations.csv` | `ea8c8312e78d1ebcea7ed5a30b5b6d91453028ad9f7d80d81a2301a8e771bd8e` |
| `WAVE_MARKERS_RM3_markers.csv` | `76b87113ae46867f988be134eb61015ef99cf544e1fe49860fe35e81cd071a97` |
| `WAVE_MARKERS_RM3_origin.csv` | `72ba87a5901a486ff473da8585d44a8b842fa53882f76f69bba4ba3325fdd1bd` |
| `WAVE_MARKERS_RM3_wave_parameters.csv` | `7308b0ae4c0cfccc45e79a381c446b75e318af29011521079d27de2d3f70cb8a` |

The archive SHA-256 is
`c767a06d2b641857dfe999cf32861a34509ba3a90fa5ced67e7f92848db90613`.
The ordinary-CI gate reproduces every marker and the origin elevation from
the published 2.5 m, 8 s, 200 m-depth sea at the unchanged `1e-10` m limit.

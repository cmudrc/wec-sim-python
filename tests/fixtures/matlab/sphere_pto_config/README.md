# Configured Sphere PTO trajectories

The `Sphere_PTO_Config` job in this [R2025b source
run](https://github.com/cmudrc/wec-sim-python/actions/runs/38472039571)
completed successfully at repo commit
`ca5064ee2126890b5c581c37e9d4b6c0ed9099cf`; unrelated jobs caused the
overall matrix run to be cancelled. The job pins WEC-Sim core
`0753b2e47f2457c078751dcfe5d251d1767b80ab` and Applications
`d53d4d4c9eda2581f04204f5d394a6ef84bb099e`.

`sphere.h5` is the exact generated source HDF5 (SHA-256
`01d11a9365aa823e971302dd4559e58d0bbf1b5e8a9841be0b61dc840f54ca7d`).
`published_400s.npz` contains all 20,001 unrounded samples of each MATLAB
body, PTO, and controller CSV for both the configured passive case and its
pretension variant, plus the source pretension parameters. Its SHA-256 is
`d4ac9704c773e49479157bb0ef085c149c30303365d29363d77b2b6d4856eb35`.
The source CSV SHA-256 hashes are:

| CSV suffix after `Sphere_PTO_Config_` | SHA-256 |
| --- | --- |
| `Configured_body1` | `35ce8e815a45e2f1e9804349acb9c50e85ac8a124cc6d973e092b0037b820e5a` |
| `Configured_pto1` | `5f78fc6de58ab0e108528985ff5a68e01ad15f83763b1db37057b21f928dce3b` |
| `controller` | `04c0f88d87a1d423249139c07eb439dcc5d6dbb16b295dd9924bd470cad24235` |
| `ConfiguredPretension_body1` | `3bba42138ef6f774d888d352aa86f26f64c03edff35a01883bf93891eb90d7d3` |
| `ConfiguredPretension_pto1` | `a1b945526c0b1e3c6d6436ce2cedf03d4d106c146bd7bb540c1391daa5702ada` |
| `ConfiguredPretension_controller` | `0ee9c1ac0a911a75e328ea7f23aed463cd90b4a33c2ef3e59470a02b1609b64d` |
| `pretension_parameters` | `0cff7f675991ab1d4abcb9a33d6ae457a5b0d5608fb5a6c297b6123ffab1dbd6` |

`test_sphere_pto_config_parity.py` uses the archive in ordinary CI and the
original CSVs when the source workflow supplies them. The independent public
`WEC` run and every existing 400 s motion, PTO, controller force, and power
limit are unchanged. These derived configurations were already documented as
passing in the parity ledger; this fixture does not add a published case badge.

# Published WaveBot impedance experiment

The [successful 600 s R2025b source run](https://github.com/cmudrc/wec-sim-python/actions/runs/38465422104)
produced `matlab-WAVEBOT_IMPEDANCE-baseline` at repo commit
`794655800211282d1f955a71f8659e11e500002c`. The workflow pins WEC-Sim
core `0753b2e47f2457c078751dcfe5d251d1767b80ab` and Applications
`d53d4d4c9eda2581f04204f5d394a6ef84bb099e`.

`waveBotBuoy.h5` is the source run's exact generated hydrodynamic file (SHA-256
`6449a44c94925dd2e149f3c0ffbd5f7941677f51989781fe762982b9f237a20a`).
`multisine3DOFA.mat` is copied without modification from the pinned
Applications `Load_Mitigating_Controls/CalcImpedance` input (SHA-256
`dba5dbc8ac72bb20054278c0ea1695c2c03784be9bb2e8af5efc0b26313888a7`).

`published_600s.npz` stores float64 values from all 60,001 published body
samples and the raw actuation time series. No source position, velocity, or
force is fed into the Python integrator. Its SHA-256 is
`df5dd29c28cf0726cfa9e0590820b2cf44a47091e31c56ab94f0229764dd72de`.
The extraction is:

| Source CSV suffix | SHA-256 | Archive data |
| --- | --- | --- |
| `CalcImpedance_body1` | `926f37cadb5e7be38738676fde63bf65a0a7bddbe8e116f9d35b3062ce91d302` | time 0, position 1/3/5, speed 7/9/11 |
| `actuation` | `858eaabffde6b63efd2453dac1d1aec74d73b5662288a2d267bc16f13538f685` | all columns |
| `mooring` | `dc9ca764ecf3fac72b949d9fcf9e9556b1eb2d1ea46d5d42c79cefd925267775` | force 13/15/17 |
| `body1_forces` | `0fd447b08f144930c1a0a8a37176294226ad5cfac9982473e144a5d3f6cc5a47` | negative damping 25/27/29 and drag 19/21/23 |

`test_wavebot_impedance_parity.py` uses this fixture in ordinary Python CI
and the full CSVs during a live MATLAB source workflow. Both paths apply the
same preexisting 600 s motion and force gates. The published input's
linearly indexed damping matrix is selected explicitly; the second test
checks that the dissipative diagonal physical default remains distinct.
This case was already counted in the WaveBot badge before the portable gate.

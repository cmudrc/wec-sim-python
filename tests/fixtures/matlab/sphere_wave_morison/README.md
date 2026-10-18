# Sphere moving Morison in regular waves, 40 s

`source.npz` retains all 4,001 float64 samples of the body's 37 logged fields,
wave elevation, and direct Morison force replay from the [passing R2025b
run](https://github.com/cmudrc/wec-sim-python/actions/runs/37773024911)
(`matlab-SPHERE_MOVING_MORISON_WAVE-baseline`, run commit
`98f397dfba1308fd7de0959f1849102b448415ea`). Its SHA-256 is
`7eb06a7ea2db78dc0eea7962497d857b2a23eb7d786163fe49faa7775a71e5a2`.
The original body, wave, and direct-force CSV SHA-256 values are
`02212b0a6670ebb03a762e926f8a0540231dcd350a4ca14f3c02514830e52a01`,
`f03e2fbc66f1a837cd345240d8a04827a0a1df41b22aab984b13b95f6c47930d`,
and `c1ceb085a6ac0e488f46cbd7d03b4a976d22eeae496bb0a2a86ae3419e313618`.

`sphere.h5` is the same run's R2025b BEMIO output, SHA-256
`e607b732510afa912dae3c7de50d848d09ae1351a008bb6ff23f3113e2b7b07c`.
The workflow pins WEC-Sim core `0753b2e47f2457c078751dcfe5d251d1767b80ab`
and Applications `d53d4d4c9eda2581f04204f5d394a6ef84bb099e`.

The ordinary Python tests independently advance the 40 s heave-only and
published surge/heave/pitch WEC configurations from the HDF5 and wave
parameters. Saved source motion and force are comparison targets, never
simulation inputs. The heave-only test retains its wave, motion, and
post-startup force limits and checks the source Morison law on projected
heave states.

The three-coordinate test now gates all 4,001 wave and active-body samples.
Maximum position differences are 12.19 mm surge, 1.003 mm heave, and
`9.78e-5` rad pitch; speed differences are 1.312/1.316 mm/s and
`4.57e-5` rad/s. After the first second, physical Morison load differences
are 429 N surge, 149 N heave, and 175 N m pitch. It uses the same explicit
limits as the earlier external-source three-coordinate test: 15/1.5 mm,
`1.5e-4` rad, 2/2 mm/s, `6e-5` rad/s, and 500/200 N plus 250 N m.
Inactive motion channels match zero. The source's startup acceleration
feedback and nonorthogonal rotation convention remain documented limits;
this strengthens an already counted Sphere case without adding a badge gate.
